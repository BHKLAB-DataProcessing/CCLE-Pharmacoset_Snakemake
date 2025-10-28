## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  if (length(snakemake@log) > 0) {
    sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
  }

  save.image(
    file.path("resources", paste0(snakemake@rule, ".RData"))
  )
} else {
  load(file.path("resources", "make_RPPA_SE.RData"))
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

rppa_cfg <- snakemake@config$molecularProfiles$rppa
rppa <- rppa_cfg$harmonized_rppa

release <- rppa_cfg$release %||% rppa$release %||% NA_character_

data_source <- rppa
data_source$release <- data_source$release %||% release
data_source$dataset <- "harmonized_rppa"

message("Starting RPPA SummarizedExperiment construction")
message("Using RPPA dataset: harmonized_rppa")

## ------------------- Load Inputs ------------------- ##
message("Reading RPPA matrix: ", INPUT$matrix)
rppa_dt <- data.table::fread(INPUT$matrix, showProgress = FALSE)

if (ncol(rppa_dt) < 2) {
  stop("RPPA matrix appears to have no feature columns.", call. = FALSE)
}

setnames(rppa_dt, 1, "depmap_id")

if (anyNA(rppa_dt$depmap_id)) {
  stop("Missing DepMap identifiers detected in RPPA matrix.", call. = FALSE)
}

message("Rows (samples): ", nrow(rppa_dt))
message("Columns (features + sample id): ", ncol(rppa_dt))

message("Reading annotated sample metadata: ", INPUT$sample_metadata)
sample_metadata <- data.table::fread(INPUT$sample_metadata, showProgress = FALSE)

if (!all(c("CCLE.depMapID", "CCLE.sampleid", "cellosaurus.cellLineName") %in% colnames(sample_metadata))) {
  stop("Annotated sample metadata is missing required columns.", call. = FALSE)
}

sample_lookup <- unique(
  sample_metadata[, .(
    CCLE.depMapID,
    CCLE.sampleid,
    cellosaurus.cellLineName,
    CCLE.site_Primary
  )],
  by = "CCLE.depMapID"
)

setkey(sample_lookup, CCLE.depMapID)

depmap_ids <- rppa_dt$depmap_id

col_meta <- sample_lookup[data.table::J(depmap_ids)]
col_meta[, depmap_id := depmap_ids]

missing_mask <- is.na(col_meta$CCLE.sampleid)
if (any(missing_mask)) {
  missing_ids <- unique(col_meta$depmap_id[missing_mask])
  message(
    "Warning: ", length(missing_ids),
    " sample(s) lack CCLE.sampleid mapping: ",
    paste(missing_ids, collapse = ", ")
  )
  # Drop samples without a mapping to keep MAE construction consistent
  keep_idx <- !missing_mask
  col_meta <- col_meta[keep_idx]
  depmap_ids <- depmap_ids[keep_idx]
  rppa_dt <- rppa_dt[keep_idx]
}

if (nrow(col_meta) == 0) {
  stop("No RPPA samples could be mapped to CCLE sample identifiers.", call. = FALSE)
}

sample_ids <- col_meta$CCLE.sampleid

# Build assay matrix (features x samples)
rppa_values <- as.matrix(rppa_dt[, -1, with = FALSE])
rownames(rppa_values) <- depmap_ids
rppa_matrix <- t(rppa_values)

if (!is.numeric(rppa_matrix[1, 1])) {
  stop("RPPA matrix failed to convert to numeric values.", call. = FALSE)
}

colnames(rppa_matrix) <- sample_ids

# Row / feature metadata -------------------------------------------------- #
feature_ids <- rownames(rppa_matrix)

feature_map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(feature_ids),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "UNIPROT"
)

feature_map_dt <- if (is.null(feature_map)) {
  data.table::data.table(UNIPROT = character())
} else {
  data.table::as.data.table(feature_map)
}

feature_summary <- feature_map_dt[
  ,
  .(
    gene_symbol = paste(sort(unique(na.omit(SYMBOL))), collapse = "|"),
    gene_name = paste(sort(unique(na.omit(GENENAME))), collapse = "|"),
    entrez_id = paste(sort(unique(na.omit(ENTREZID))), collapse = "|"),
    n_mapped_genes = uniqueN(na.omit(SYMBOL))
  ),
  by = UNIPROT
]

if (nrow(feature_summary)) {
  feature_summary[
    gene_symbol == "",
    gene_symbol := NA_character_
  ]
  feature_summary[
    gene_name == "",
    gene_name := NA_character_
  ]
  feature_summary[
    entrez_id == "",
    entrez_id := NA_character_
  ]
}

feature_metadata <- data.table(uniprot_id = feature_ids)[
  feature_summary,
  on = c("uniprot_id" = "UNIPROT")
]

feature_metadata[
  is.na(n_mapped_genes),
  n_mapped_genes := 0L
]

setcolorder(feature_metadata, c("uniprot_id", "gene_symbol", "gene_name", "entrez_id", "n_mapped_genes"))

# Column metadata --------------------------------------------------------- #

col_metadata <- data.table(
  sampleid = sample_ids,
  depmap_id = col_meta$CCLE.depMapID,
  cell_line = col_meta$cellosaurus.cellLineName,
  tissue = col_meta$CCLE.site_Primary
)

col_metadata[
  is.na(cell_line),
  cell_line := "unknown"
]

col_metadata[
  is.na(tissue),
  tissue := "unknown"
]

# SummarizedExperiment ---------------------------------------------------- #

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = rppa_matrix),
  rowData = S4Vectors::DataFrame(feature_metadata, row.names = feature_metadata$uniprot_id),
  colData = S4Vectors::DataFrame(
    sampleid = col_metadata$sampleid,
    depmap_id = col_metadata$depmap_id,
    cell_line = col_metadata$cell_line,
    tissue = col_metadata$tissue,
    batchid = rep(NA, nrow(col_metadata)),
    row.names = col_metadata$sampleid
  )
)

SummarizedExperiment::metadata(se) <- list(
  annotation = "proteomics",
  datatype = "proteomics",
  class = "SummarizedExperiment",
  filename = basename(INPUT$matrix),
  release = release,
  data_source = data_source,
  numSamples = ncol(se),
  numFeatures = nrow(se),
  generated_on = as.character(Sys.Date())
)

message("RPPA SummarizedExperiment summary:")
show(se)

# Outputs ----------------------------------------------------------------- #

dir.create(dirname(OUTPUT$rppa_se), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT$feature_metadata), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT$sample_metadata), recursive = TRUE, showWarnings = FALSE)

se_list <- list(
  proteomics.rppa = se
)

saveRDS(se_list, file = OUTPUT$rppa_se)
data.table::fwrite(feature_metadata, file = OUTPUT$feature_metadata, sep = "\t")
data.table::fwrite(col_metadata, file = OUTPUT$sample_metadata, sep = "\t")

message("Saved outputs:")
message(" - SummarizedExperiment list: ", OUTPUT$rppa_se)
message(" - Feature metadata: ", OUTPUT$feature_metadata)
message(" - Sample metadata: ", OUTPUT$sample_metadata)

if (exists("snakemake")) {
  sink()
}
