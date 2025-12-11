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
  if (file.exists(file.path("snapshots/", paste0(snakemake@rule, ".RData")))) {
    load(file.path("snapshots/", paste0(snakemake@rule, ".RData")))
  }
}
mutation_cache <- file.path("resources/make_Mutation_SE.RData")
if (file.exists(mutation_cache)) {
  load(mutation_cache)
}

message("Loading mutation data")
input <- readRDS(INPUT$preprocessedMutation)
assay <- input$assay
gene_meta <- data.table::as.data.table(input$gene_meta)

message("Loading GENCODE annotation for ranges")
gencode_gr <- rtracklayer::import(INPUT$GENCODE_Annotation)
gencode_genes <- gencode_gr[gencode_gr$type == "gene", ]
gencode_dt <- data.table::data.table(
  gene_symbol = gencode_genes$gene_name,
  gene_id = sub("\\.\\d+$", "", gencode_genes$gene_id),
  seqnames = as.character(GenomicRanges::seqnames(gencode_genes)),
  start = BiocGenerics::start(gencode_genes),
  end = BiocGenerics::end(gencode_genes)
)
data.table::setkey(gencode_dt, gene_id)
gencode_by_symbol <- data.table::copy(gencode_dt)
data.table::setkey(gencode_by_symbol, gene_symbol)

lookup_ids <- sub("\\.\\d+$", "", gene_meta$EnsemblGeneID)
gene_meta[, EnsemblGeneID := lookup_ids]

match_by_id <- gencode_dt[lookup_ids]
missing_id <- which(is.na(match_by_id$gene_id))

if (length(missing_id)) {
  sym_hits <- gencode_by_symbol[gene_meta$HugoSymbol[missing_id]]
  replaceable <- which(!is.na(sym_hits$gene_id))
  match_by_id[missing_id[replaceable]] <- sym_hits[replaceable]
}

available <- !is.na(match_by_id$gene_id)
if (!all(available)) {
  dropped <- gene_meta$HugoSymbol[!available]
  preview <- head(dropped, 20)
  warning(sprintf(
    "Dropping %d genes without GENCODE mapping (showing up to 20): %s",
    length(dropped),
    paste(preview, collapse = ", ")
  ))
}

match_by_id <- match_by_id[available]
assay <- assay[available, , drop = FALSE]
gene_meta <- gene_meta[available]

complete_coords <- !is.na(match_by_id$start) & !is.na(match_by_id$end)
if (!all(complete_coords)) {
  dropped <- gene_meta$HugoSymbol[!complete_coords]
  preview <- head(dropped, 20)
  warning(sprintf(
    "Dropping %d genes with missing coordinates (showing up to 20): %s",
    length(dropped),
    paste(preview, collapse = ", ")
  ))
}

match_by_id <- match_by_id[complete_coords]
assay <- assay[complete_coords, , drop = FALSE]
gene_meta <- gene_meta[complete_coords]

message("Creating gene-level rowRanges from GENCODE")
rowRanges <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(match_by_id$seqnames),
  ranges = IRanges::IRanges(start = match_by_id$start, end = match_by_id$end),
  strand = S4Vectors::Rle("*"),
  gene_symbol = gene_meta$HugoSymbol,
  EntrezGeneID = gene_meta$EntrezGeneID,
  EnsemblGeneID = gene_meta$EnsemblGeneID
)
names(rowRanges) <- gene_meta$HugoSymbol

message("Creating mutation SummarizedExperiment")
mutation_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    exprs = assay
  ),
  rowRanges = rowRanges,
  colData = data.table::data.table(
    sampleid = colnames(assay),
    batchid = rep(NA, ncol(assay))
  )
)

message("Adding metadata to mutation SummarizedExperiment")
(mutation_se@metadata <- list(
  annotation = "mut",
  datatype = "genes",
  class = "RangedSummarizedExperiment",
  filename = basename(INPUT$preprocessedMutation),
  data_source = snakemake@config$molecularProfiles$mutation$somatic_mutations,
  date = Sys.Date(),
  numSamples = ncol(mutation_se),
  numGenes = nrow(mutation_se)
))

show(mutation_se)

se_list <- list(
  mut.genes = mutation_se
)

message("Saving mutation SE to ", OUTPUT$processedMutationSE)
saveRDS(se_list, file = OUTPUT$processedMutationSE)
