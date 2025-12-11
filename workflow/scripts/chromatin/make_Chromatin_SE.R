## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  PARAMS <- snakemake@params
  LOG <- snakemake@log

  if (length(LOG)) {
    sink(LOG[[1]], split = TRUE, type = c("output", "message"))
  }

  dir.create("snapshots", showWarnings = FALSE, recursive = TRUE)
  save.image(file.path("snapshots", paste0(snakemake@rule, ".RData")))
} else {
  if (file.exists(file.path("snapshots", "make_Chromatin_SE.RData"))) {
    load(file.path("snapshots", "make_Chromatin_SE.RData"))
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(S4Vectors)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

chrom_cfg <- PARAMS$dataset_config
dataset_label <- "chromatin.global_mrm"

message("Building chromatin (MRM) SummarizedExperiment")
message(" - Input matrix: ", INPUT$data)
message(" - Annotated sample metadata: ", INPUT$sample_metadata)

dt <- data.table::fread(INPUT$data, showProgress = FALSE)

if (!all(c("CellLineName", "BroadID") %in% names(dt))) {
  stop("Chromatin file is missing CellLineName and BroadID columns.", call. = FALSE)
}

feature_cols <- setdiff(names(dt), c("CellLineName", "BroadID"))
if (!length(feature_cols)) {
  stop("No feature columns found in chromatin file.", call. = FALSE)
}

# Sample metadata join ---------------------------------------------------- #
sample_md <- data.table::fread(INPUT$sample_metadata, showProgress = FALSE)

required_cols <- c("CCLE.depMapID", "CCLE.sampleid", "cellosaurus.cellLineName")
if (!all(required_cols %in% names(sample_md))) {
  stop("Annotated sample metadata missing required columns.", call. = FALSE)
}

sample_lookup <- unique(sample_md[, .(
  CCLE.depMapID,
  CCLE.sampleid,
  cellosaurus.cellLineName,
  CCLE.site_Primary,
  CCLE.type
)])

col_meta <- data.table::data.table(
  BroadID = dt$BroadID,
  dataset_sample_id = dt$CellLineName
)

setkey(sample_lookup, CCLE.depMapID)
col_meta <- sample_lookup[col_meta, on = c("CCLE.depMapID" = "BroadID")]

missing_mask <- is.na(col_meta$CCLE.sampleid)
if (any(missing_mask)) {
  setkey(sample_lookup, cellosaurus.cellLineName)
  col_meta[missing_mask, c(
    "CCLE.depMapID",
    "CCLE.sampleid",
    "cellosaurus.cellLineName",
    "CCLE.site_Primary",
    "CCLE.type"
  ) := sample_lookup[
    dataset_sample_id,
    .(
      CCLE.depMapID,
      CCLE.sampleid,
      cellosaurus.cellLineName,
      CCLE.site_Primary,
      CCLE.type
    ),
    on = "cellosaurus.cellLineName"
  ]
  ]
}

if (any(is.na(col_meta$CCLE.sampleid))) {
  unresolved <- unique(col_meta$dataset_sample_id[is.na(col_meta$CCLE.sampleid)])
  stop("Failed to map chromatin samples to CCLE.sampleid: ",
       paste(unresolved, collapse = ", "))
}

col_meta[, dataset_sample_id := dt$CellLineName]
col_meta[, BroadID := dt$BroadID]

# Assay matrix ------------------------------------------------------------ #
assay_mat <- as.matrix(dt[, ..feature_cols])
storage.mode(assay_mat) <- "numeric"
rownames(assay_mat) <- dt$CCLE.sampleid <- col_meta$CCLE.sampleid
assay_mat <- t(assay_mat)  # features x samples

stopifnot(!any(is.na(colnames(assay_mat))))

# Feature metadata
feature_md <- data.table::data.table(
  feature_id = feature_cols,
  histone_mark = feature_cols
)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(intensity = assay_mat),
  colData = S4Vectors::DataFrame(col_meta[, .(
    sampleid = CCLE.sampleid,
    dataset_sample_id,
    depmap_id = CCLE.depMapID,
    broad_id = BroadID,
    tissue = CCLE.site_Primary,
    type = CCLE.type
  )], row.names = col_meta$CCLE.sampleid),
  rowData = S4Vectors::DataFrame(feature_md, row.names = feature_md$feature_id),
  metadata = list(
    annotation = "chromatin",
    datatype = "global_histone_mrm",
    class = "SummarizedExperiment",
    data_source = chrom_cfg,
    date = chrom_cfg$release %||% Sys.Date(),
    numSamples = ncol(assay_mat),
    numFeatures = nrow(assay_mat)
  )
)

message("Samples mapped: ", ncol(se))
message("Features: ", nrow(se))

# Save outputs ------------------------------------------------------------ #
dir.create(dirname(OUTPUT$chromatin_se), showWarnings = FALSE, recursive = TRUE)

saveRDS(se, file = OUTPUT$chromatin_se)
data.table::fwrite(feature_md, OUTPUT$feature_metadata, sep = "\t")
data.table::fwrite(as.data.table(colData(se)), OUTPUT$sample_metadata, sep = "\t")

message("Chromatin SummarizedExperiment written to ", OUTPUT$chromatin_se)

if (length(LOG)) sink()
