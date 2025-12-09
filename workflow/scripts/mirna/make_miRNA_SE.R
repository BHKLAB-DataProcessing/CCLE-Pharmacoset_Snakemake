## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  LOG <- snakemake@log

  if (length(LOG)) {
    sink(LOG[[1]], split = TRUE, type = c("output", "message"))
  }

  dir.create("snapshots", showWarnings = FALSE, recursive = TRUE)
  save.image(file.path("snapshots", paste0(snakemake@rule, ".RData")))
} else {
  if (file.exists(file.path("snapshots", "make_miRNA_SE.RData"))) {
    load(file.path("snapshots", "make_miRNA_SE.RData"))
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
})

read_gct <- function(path) {
  # GCT format has two header rows: version and dims. Skip them.
  dt <- data.table::fread(path, skip = 2)
  value_cols <- setdiff(colnames(dt), c("Name", "Description"))
  mat <- as.matrix(dt[, ..value_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- dt[["Name"]]
  rowData <- DataFrame(feature_id = dt[["Name"]], description = dt[["Description"]])
  colData <- DataFrame(sampleid = value_cols, row.names = value_cols)
  SummarizedExperiment(
    assays = list(exprs = mat),
    rowData = rowData,
    colData = colData,
    metadata = list(annotation = "mirna", datatype = "mirna_gct")
  )
}

message("Reading miRNA GCT: ", INPUT$gct)
mirna_gct_se <- read_gct(INPUT$gct)

message("Reading miRNA MIMAT CSV: ", INPUT$mimat)
mimat_dt <- data.table::fread(INPUT$mimat)
stopifnot(ncol(mimat_dt) > 1)
# Format: rows = samples, columns = MIMAT IDs, first column blank/sampleid
sample_ids <- mimat_dt[[1]]
feature_ids <- colnames(mimat_dt)[-1]
mat <- as.matrix(mimat_dt[, -1, with = FALSE])
storage.mode(mat) <- "numeric"
mat <- t(mat) # features x samples
rownames(mat) <- feature_ids
colnames(mat) <- sample_ids
rowData <- DataFrame(feature_id = feature_ids)
colData <- DataFrame(sampleid = sample_ids, row.names = sample_ids)
mirna_mimat_se <- SummarizedExperiment(
  assays = list(exprs = mat),
  rowData = rowData,
  colData = colData,
  metadata = list(annotation = "mirna", datatype = "mirna_mimat")
)

message("Saving outputs")
saveRDS(mirna_gct_se, file = OUTPUT$mirna_gct_se)
saveRDS(mirna_mimat_se, file = OUTPUT$mirna_mimat_se)

if (length(LOG)) sink()
