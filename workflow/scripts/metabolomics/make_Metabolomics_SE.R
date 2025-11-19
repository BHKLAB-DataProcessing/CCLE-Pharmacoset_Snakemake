## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  PARAMS <- snakemake@params
  THREADS <- snakemake@threads

  if (length(snakemake@log)) {
    sink(
      snakemake@log[[1]],
      append = FALSE,
      type = c("output", "message"),
      split = TRUE
    )
  }

  dir.create("snapshots", showWarnings = FALSE, recursive = TRUE)
  save.image(
    file.path("snapshots/", paste0(snakemake@rule, ".RData"))
  )
} else {
  if (file.exists("snapshots/make_Metabolomics_SE.RData")) {
    load("snapshots/make_Metabolomics_SE.RData")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(S4Vectors)
})

## ------------------- Helper Functions ------------------- ##
clean_feature_id <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- "metabolite"
  make.unique(x)
}

## ------------------- Load Inputs ------------------- ##
message("Reading metabolomics matrix: ", INPUT$data)
met_dt <- data.table::fread(
  INPUT$data,
  na.strings = c("NA", "NaN", ""),
  showProgress = FALSE
)

if (!nrow(met_dt)) {
  stop("Metabolomics table is empty.", call. = FALSE)
}

required_cols <- c("CCLE_ID", "DepMap_ID")
missing_cols <- setdiff(required_cols, names(met_dt))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

expr_cols <- setdiff(names(met_dt), required_cols)
if (!length(expr_cols)) {
  stop("No metabolite columns detected.", call. = FALSE)
}

sample_ids <- met_dt$CCLE_ID
depmap_ids <- met_dt$DepMap_ID

if (anyDuplicated(sample_ids)) {
  dupes <- unique(sample_ids[duplicated(sample_ids)])
  stop(
    "Duplicate CCLE_ID entries detected: ",
    paste(head(dupes, 10), collapse = ", ")
  )
}

message("Metabolites detected: ", length(expr_cols))
message("Samples detected: ", length(sample_ids))

expr_matrix <- t(as.matrix(met_dt[, ..expr_cols]))
storage.mode(expr_matrix) <- "numeric"
colnames(expr_matrix) <- sample_ids

feature_ids <- clean_feature_id(expr_cols)
rownames(expr_matrix) <- feature_ids

feature_dt <- data.table::data.table(
  feature_id = feature_ids,
  metabolite = expr_cols,
  measurement = "relative_abundance",
  platform = "LC-MS"
)

## ------------------- Sample Metadata Alignment ------------------- ##
sample_metadata <- data.table::fread(
  INPUT$sample_metadata,
  na.strings = "NA",
  sep = "\t",
  showProgress = FALSE
)

sample_metadata <- sample_metadata[!duplicated(CCLE.sampleid)]
sample_metadata[, sampleid := CCLE.sampleid]
setkey(sample_metadata, sampleid)

column_map <- data.table::data.table(
  dataset_sample_id = sample_ids,
  depmap_dataset_id = depmap_ids,
  sampleid = sample_ids
)

column_map <- merge(
  column_map,
  sample_metadata[,
    .(
      sampleid,
      depmap_id = CCLE.depMapID,
      cell_line = CCLE.name,
      tissue = CCLE.site_Primary
    )
  ],
  by = "sampleid",
  all.x = TRUE,
  sort = FALSE
)

unmapped <- column_map[is.na(depmap_id)]
if (nrow(unmapped)) {
  message(
    sprintf(
      "Dropping %d samples without annotated metadata. Examples: %s",
      nrow(unmapped),
      paste(head(unmapped$dataset_sample_id, 10), collapse = ", ")
    )
  )
}

column_map <- column_map[!is.na(depmap_id)]

if (!nrow(column_map)) {
  stop("No metabolomics samples matched the annotated metadata.", call. = FALSE)
}

keep_cols <- column_map$dataset_sample_id
expr_matrix <- expr_matrix[, keep_cols, drop = FALSE]

col_metadata_dt <- data.table::data.table(
  sampleid = column_map$sampleid,
  depmap_id = column_map$depmap_id,
  dataset_sample_id = column_map$dataset_sample_id,
  depmap_dataset_id = column_map$depmap_dataset_id,
  cell_line = column_map$cell_line,
  tissue = column_map$tissue,
  batchid = NA_character_
)

colData <- S4Vectors::DataFrame(
  sampleid = col_metadata_dt$sampleid,
  depmap_id = col_metadata_dt$depmap_id,
  dataset_sample_id = col_metadata_dt$dataset_sample_id,
  depmap_dataset_id = col_metadata_dt$depmap_dataset_id,
  cell_line = col_metadata_dt$cell_line,
  tissue = col_metadata_dt$tissue,
  batchid = col_metadata_dt$batchid,
  row.names = col_metadata_dt$sampleid
)

rowData <- S4Vectors::DataFrame(
  feature_dt,
  row.names = feature_dt$feature_id
)

metabolomics_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = expr_matrix),
  rowData = rowData,
  colData = colData
)

SummarizedExperiment::metadata(metabolomics_se) <- c(
  list(
    annotation = "metabolomics",
    datatype = "metabolites",
    class = "SummarizedExperiment",
    filename = basename(INPUT$data),
    data_source = PARAMS$dataset_config,
    numSamples = ncol(metabolomics_se),
    numFeatures = nrow(metabolomics_se),
    generated_on = as.character(Sys.Date())
  )
)

message("Metabolomics SummarizedExperiment summary:")
show(metabolomics_se)

## ------------------- Write Outputs ------------------- ##
dir.create(
  dirname(OUTPUT$metabolomics_se),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(metabolomics.lcms = metabolomics_se),
  file = OUTPUT$metabolomics_se
)

data.table::fwrite(
  feature_dt,
  file = OUTPUT$feature_metadata,
  sep = "\t"
)

data.table::fwrite(
  col_metadata_dt,
  file = OUTPUT$sample_metadata,
  sep = "\t"
)

message("Saved metabolomics SE and metadata:")
message(" - SummarizedExperiment list: ", OUTPUT$metabolomics_se)
message(" - Feature metadata: ", OUTPUT$feature_metadata)
message(" - Sample metadata: ", OUTPUT$sample_metadata)

sink()
