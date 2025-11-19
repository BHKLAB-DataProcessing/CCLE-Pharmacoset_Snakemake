## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  # setup logger if log file is provided
  if (length(snakemake@log) > 0) {
    sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
  }

  dir.create("snapshots", showWarnings = FALSE, recursive = TRUE)
  save.image(
    file.path("snapshots/", paste0(snakemake@rule, ".RData"))
  )
} else {
  if (file.exists(file.path("snapshots/make_CNV_SE.RData"))) {
    load(file.path("snapshots/make_CNV_SE.RData"))
  }
}

library(GenomicRanges)
library(SummarizedExperiment)

## ------------------- LOAD INPUT ------------------- ##

message("Loading CNV data")
cnvDt <- data.table::fread(
  INPUT$cnv,
  header = TRUE,
  stringsAsFactors = FALSE
)

if ("Unnamed: 0" %in% names(cnvDt)) {
  cnvDt[, "Unnamed: 0" := NULL]
}

metadata_cols <- intersect(
  c(
    "SequencingID",
    "ModelID",
    "IsDefaultEntryForModel",
    "ModelConditionID",
    "IsDefaultEntryForMC"
  ),
  names(cnvDt)
)

message("Filtering to default model entries")
cnvDt[, IsDefaultEntryForModel := tolower(IsDefaultEntryForModel)]
cnvDt[, IsDefaultEntryForMC := tolower(IsDefaultEntryForMC)]
cnv_filtered <- cnvDt[
  IsDefaultEntryForModel == "yes" &
    (is.na(IsDefaultEntryForMC) | IsDefaultEntryForMC == "yes")
]
if (!nrow(cnv_filtered)) {
  warning(
    "No rows marked as default for both model and condition; using all rows"
  )
  cnv_filtered <- cnvDt
}

dupe_models <- cnv_filtered[duplicated(ModelID), unique(ModelID)]
if (length(dupe_models)) {
  warning(
    "Duplicate ModelID entries detected, keeping first occurrence for: ",
    paste(dupe_models, collapse = ", ")
  )
  cnv_filtered <- cnv_filtered[!duplicated(ModelID)]
}

gene_cols <- setdiff(names(cnv_filtered), metadata_cols)

message("Preparing gene annotations from column names")
extract_ensembl <- function(x) {
  match <- regexpr("ENSG[0-9]+(?:\\.[0-9]+)?", x)
  ifelse(match > 0, regmatches(x, match), NA_character_)
}

gene_annotation <- data.table::data.table(
  original = gene_cols,
  gene_symbol = trimws(sub("^\\s*([^\\s(]+).*", "\\1", gene_cols)),
  ensembl_id = vapply(gene_cols, extract_ensembl, character(1))
)
gene_annotation[, ensembl_id := sub("\\.\\d+$", "", ensembl_id)]
gene_annotation <- unique(gene_annotation, by = "original")

message("Loading sample metadata for ModelID to CCLE.sampleid mapping")
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)
sample_lookup <- unique(sampleMetadata[!is.na(CCLE.depMapID), .(
  ModelID = CCLE.depMapID,
  CCLE.sampleid
)])
data.table::setkey(sample_lookup, ModelID)

sample_matches <- sample_lookup[cnv_filtered$ModelID]
missing_samples <- sample_matches[is.na(CCLE.sampleid), ModelID]
if (length(missing_samples)) {
  preview <- head(missing_samples, 20)
  warning(
    sprintf(
      "Using ModelID as sample name for %d samples with no CCLE.sampleid match (showing up to 20): %s",
      length(missing_samples),
      paste(preview, collapse = ", ")
    )
  )
}
sample_matches[, sampleid_final := data.table::fifelse(
  is.na(CCLE.sampleid),
  ModelID,
  CCLE.sampleid
)]

message("Building numeric CNV matrix")
cnv_matrix <- as.matrix(cnv_filtered[, ..gene_cols])
storage.mode(cnv_matrix) <- "numeric"
colnames(cnv_matrix) <- gene_annotation$gene_symbol
rownames(cnv_matrix) <- sample_matches$ModelID

if (anyDuplicated(colnames(cnv_matrix))) {
  message("Collapsing duplicated gene symbols by mean")
  unique_syms <- unique(colnames(cnv_matrix))
  collapsed <- lapply(unique_syms, function(sym) {
    cols <- which(colnames(cnv_matrix) == sym)
    if (length(cols) == 1) cnv_matrix[, cols]
    else rowMeans(cnv_matrix[, cols, drop = FALSE], na.rm = TRUE)
  })
  cnv_matrix <- do.call(cbind, collapsed)
  colnames(cnv_matrix) <- unique_syms
}

cnv_matrix <- t(cnv_matrix)
colnames(cnv_matrix) <- sample_matches$sampleid_final

message("Loading GENCODE annotation")
gencode_gr <- rtracklayer::import(INPUT$GENCODE_Annotation)
gencode_genes <- gencode_gr[gencode_gr$type == "gene", ]
gencode_dt <- data.table::data.table(
  gene_symbol = gencode_genes$gene_name,
  gene_id = sub("\\.\\d+$", "", gencode_genes$gene_id),
  seqnames = as.character(GenomicRanges::seqnames(gencode_genes)),
  start = BiocGenerics::start(gencode_genes),
  end = BiocGenerics::end(gencode_genes)
)
gencode_dt <- unique(gencode_dt, by = "gene_symbol")
data.table::setkey(gencode_dt, gene_symbol)

gene_symbols <- rownames(cnv_matrix)
gene_symbols <- trimws(gene_symbols)
match_idx <- match(gene_symbols, gencode_dt$gene_symbol)
matched <- which(!is.na(match_idx))

if (!length(matched)) {
  stop("No gene symbols from CNV file matched GENCODE gene names; check annotation.")
}

if (length(matched) < length(gene_symbols)) {
  dropped_genes <- gene_symbols[-matched]
  preview <- head(dropped_genes, 20)
  warning(
    sprintf(
      "Dropping %d genes without GENCODE mapping (showing up to 20): %s",
      length(dropped_genes),
      paste(preview, collapse = ", ")
    )
  )
}

symbol_match <- gencode_dt[match_idx[matched]]
symbol_match[, gene_symbol := gene_symbols[matched]]
cnv_matrix <- cnv_matrix[symbol_match$gene_symbol, , drop = FALSE]

message("Creating GRanges object for CNV data")
cnv_gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(symbol_match$seqnames),
  ranges = IRanges::IRanges(start = symbol_match$start, end = symbol_match$end),
  strand = S4Vectors::Rle("*"),
  gene_id = symbol_match$gene_id,
  gene_name = symbol_match$gene_symbol
)
names(cnv_gr) <- symbol_match$gene_symbol

message("Creating SummarizedExperiment object for CNV data")
cnv_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    exprs = cnv_matrix
  ),
  rowRanges = cnv_gr,
  colData = data.table::data.table(
    sampleid = colnames(cnv_matrix),
    dataset_sample_id = sample_matches$ModelID,
    batchid = rep(NA, ncol(cnv_matrix))
  )
)

message("Adding metadata to CNV SummarizedExperiment object")
(cnv_se@metadata <- list(
  annotation = "cnv",
  datatype = "genes",
  class = "RangedSummarizedExperiment",
  filename = basename(INPUT$cnv),
  data_source = snakemake@config$molecularProfiles$cnv$copynumber_gene_wgs,
  date = Sys.Date(),
  numSamples = ncol(cnv_se),
  numGenes = nrow(cnv_se)
))

print(cnv_se)

## ------------------- SAVE OUTPUT ------------------- ##
message("Saving CNV SummarizedExperiment object")
se_list <- list(
  cnv.genes = cnv_se
)

saveRDS(se_list, file = OUTPUT$CNV_SE)
