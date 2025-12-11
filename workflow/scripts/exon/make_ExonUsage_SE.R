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
  if (file.exists(file.path("snapshots", "make_ExonUsage_SE.RData"))) {
    load(file.path("snapshots", "make_ExonUsage_SE.RData"))
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
})

read_gct_matrix <- function(path) {
  dt <- data.table::fread(path, skip = 2)
  name_col <- colnames(dt)[1]
  desc_col <- if ("Description" %in% colnames(dt)) "Description" else NULL
  value_cols <- setdiff(colnames(dt), c("NAME", "Name", desc_col))
  mat <- as.matrix(dt[, ..value_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- dt[[name_col]]
  rowData <- DataFrame(feature_id = dt[[name_col]])
  if (!is.null(desc_col)) rowData$description <- dt[[desc_col]]
  list(mat = mat, rowData = rowData)
}

message("Reading exon usage ratio: ", INPUT$ratio)
ratio <- read_gct_matrix(INPUT$ratio)
message("Reading exon denominator counts: ", INPUT$denom)
denom <- read_gct_matrix(INPUT$denom)

colData_ratio <- DataFrame(sampleid = colnames(ratio$mat), row.names = colnames(ratio$mat))
colData_denom <- DataFrame(sampleid = colnames(denom$mat), row.names = colnames(denom$mat))

exon_ratio_se <- SummarizedExperiment(
  assays = list(ratio = ratio$mat),
  rowData = ratio$rowData,
  colData = colData_ratio,
  metadata = list(annotation = "exon", datatype = "exon_usage_ratio")
)

exon_denom_se <- SummarizedExperiment(
  assays = list(counts = denom$mat),
  rowData = denom$rowData,
  colData = colData_denom,
  metadata = list(annotation = "exon", datatype = "exon_usage_denominator")
)

saveRDS(exon_ratio_se, file = OUTPUT$exon_ratio_se)
saveRDS(exon_denom_se, file = OUTPUT$exon_denom_se)

if (length(LOG)) sink()
