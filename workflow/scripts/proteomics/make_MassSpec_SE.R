## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  PARAMS <- snakemake@params
  LOG <- snakemake@log

  if (length(LOG)) {
    sink(LOG[[1]], split = TRUE, type = c("output", "message"))
  }

  save.image(
    file.path("resources", paste0(snakemake@rule, ".RData"))
  )
} else {
  load(file.path("resources", "make_MassSpec_SE.RData"))
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(readxl)
})

`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

first_non_missing <- function(x) {
  vals <- unique(na.omit(x))
  if (length(vals)) {
    vals[1]
  } else {
    NA
  }
}

cfg <- PARAMS$dataset_config
release <- PARAMS$dataset_release %||% cfg$release %||% NA_character_

message("Building CCLE mass spectrometry SummarizedExperiment")
message(" - Quantitation file: ", INPUT$quant)
message(" - Sample information file: ", INPUT$sample_info)
message(" - Annotated sample metadata: ", INPUT$sample_metadata)

## ------------------- Load Inputs ------------------- ##

message("Reading mass spectrometry quantitation matrix")
quant_dt <- data.table::fread(
  INPUT$quant,
  sep = ",",
  showProgress = FALSE
)

if (!nrow(quant_dt)) {
  stop("Mass spectrometry quantitation table is empty.", call. = FALSE)
}

message("Rows (features): ", nrow(quant_dt))
message("Columns (all): ", ncol(quant_dt))

expr_col_mask <- grepl("_TenPx[0-9]+$", names(quant_dt))

if (!any(expr_col_mask)) {
  stop(
    "Unable to identify expression columns (expected *_TenPx###).",
    call. = FALSE
  )
}

expr_cols <- names(quant_dt)[expr_col_mask]
metadata_cols <- setdiff(names(quant_dt), expr_cols)

message("Detected ", length(expr_cols), " assay columns.")

feature_dt <- data.table::copy(quant_dt[, ..metadata_cols])

# Build feature identifiers ------------------------------------------------ #
feature_ids <- if ("Protein_Id" %in% names(feature_dt)) {
  as.character(feature_dt$Protein_Id)
} else {
  rep(NA_character_, nrow(feature_dt))
}

if ("Uniprot_Acc" %in% names(feature_dt)) {
  missing_idx <- is.na(feature_ids) | feature_ids == ""
  feature_ids[missing_idx] <- as.character(feature_dt$Uniprot_Acc[missing_idx])
}

if ("Gene_Symbol" %in% names(feature_dt)) {
  missing_idx <- is.na(feature_ids) | feature_ids == ""
  feature_ids[missing_idx] <- as.character(feature_dt$Gene_Symbol[missing_idx])
}

if (any(is.na(feature_ids) | feature_ids == "")) {
  fallback_ids <- sprintf("feature_%05d", seq_len(nrow(feature_dt)))
  missing_idx <- is.na(feature_ids) | feature_ids == ""
  feature_ids[missing_idx] <- fallback_ids[missing_idx]
}

feature_ids <- make.unique(feature_ids)

feature_dt[, feature_id := feature_ids]
setcolorder(
  feature_dt,
  c("feature_id", setdiff(names(feature_dt), "feature_id"))
)

expr_matrix <- as.matrix(quant_dt[, ..expr_cols])
storage.mode(expr_matrix) <- "numeric"
rownames(expr_matrix) <- feature_ids

rm(quant_dt)  # free memory

# Column metadata ---------------------------------------------------------- #

message("Parsing dataset column metadata")
column_map <- data.table::data.table(
  dataset_sample_id = expr_cols,
  column_index = seq_along(expr_cols)
)

column_map$sample_base <- sub(
  "_TenPx[0-9]+$",
  "",
  column_map$dataset_sample_id
)

column_map$tenplex_id <- sub(
  "^.*_TenPx([0-9]+)$",
  "\\1",
  column_map$dataset_sample_id
)

invalid_pattern <- !grepl(
  "_TenPx[0-9]+$",
  column_map$dataset_sample_id
)
if (any(invalid_pattern)) {
  column_map$sample_base[invalid_pattern] <- column_map$dataset_sample_id[invalid_pattern]
  column_map$tenplex_id[invalid_pattern] <- NA_character_
}

column_map$`Protein.10.Plex.ID` <- column_map$tenplex_id

message("Reading Gygi lab sample information workbook")
sample_info <- tryCatch(
  readxl::read_excel(
    INPUT$sample_info,
    sheet = "Sample_Information"
  ),
  error = function(e) {
    message("Failed to read 'Sample_Information' sheet: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(sample_info)) {
  sample_info_dt <- data.table::as.data.table(sample_info)
  if ("Protein.10.Plex.ID" %in% names(sample_info_dt)) {
    sample_info_dt[
      ,
      `Protein.10.Plex.ID` := as.character(`Protein.10.Plex.ID`)
    ]
  }

  old_cols <- intersect(
    names(sample_info_dt),
    c("CCLE.code", "ccle.code")
  )
  if (length(old_cols)) {
    data.table::setnames(
      sample_info_dt,
      old = old_cols,
      new = rep("CCLE.Code", length(old_cols))
    )
  }

  if ("CCLE.Code" %in% names(sample_info_dt)) {
    column_map[, CCLE.Code := sample_base]

    suppressWarnings({
      column_map <- merge(
        column_map,
        sample_info_dt,
        by = c("CCLE.Code", "Protein.10.Plex.ID"),
        all.x = TRUE,
        sort = FALSE,
        suffixes = c("", ".info")
      )
    })
  } else {
    message("Sample information file does not contain 'CCLE.Code'; proceeding without merge.")
  }
} else {
  message("Proceeding without supplemental sample information.")
}

# Annotated sample metadata ------------------------------------------------ #

message("Loading annotated sample metadata")
annotated_sample_metadata <- data.table::fread(
  INPUT$sample_metadata,
  sep = "\t",
  showProgress = FALSE
)

sample_lookup <- unique(annotated_sample_metadata[
  ,
  .(
    CCLE.sampleid,
    CCLE.depMapID,
    cellosaurus.cellLineName,
    CCLE.site_Primary
  )
])

setnames(
  sample_lookup,
  old = c("CCLE.sampleid", "CCLE.depMapID", "cellosaurus.cellLineName", "CCLE.site_Primary"),
  new = c("sampleid", "depmap_id", "cell_line", "tissue")
)

column_map[
  ,
  sampleid := sample_base
]

column_map <- merge(
  column_map,
  sample_lookup,
  by = "sampleid",
  all.x = TRUE,
  sort = FALSE
)

setorder(column_map, column_index)

unmapped <- column_map[is.na(depmap_id)]

if (nrow(unmapped)) {
  message(
    sprintf(
      "Dropping %d columns that do not map to annotated sample metadata.",
      nrow(unmapped)
    )
  )
  message(
    "Unmapped dataset identifiers: ",
    paste(head(unmapped$dataset_sample_id, 10), collapse = ", "),
    if (nrow(unmapped) > 10) " ..."
  )
}

column_map <- column_map[!is.na(depmap_id)]

if (!nrow(column_map)) {
  stop(
    "No assay columns mapped to annotated sample metadata.",
    call. = FALSE
  )
}

expr_matrix <- expr_matrix[, column_map$dataset_sample_id, drop = FALSE]

# Aggregate technical replicates ------------------------------------------ #

message("Aggregating technical replicates by CCLE sample identifier")
unique_samples <- unique(column_map$sampleid)

aggregated_matrix <- matrix(
  NA_real_,
  nrow = nrow(expr_matrix),
  ncol = length(unique_samples),
  dimnames = list(rownames(expr_matrix), unique_samples)
)

col_metadata_list <- vector("list", length(unique_samples))

for (idx in seq_along(unique_samples)) {
  sample_name <- unique_samples[[idx]]
  sample_rows <- column_map[sampleid == sample_name]
  sample_cols <- sample_rows$dataset_sample_id
  sample_matrix <- expr_matrix[, sample_cols, drop = FALSE]

  if (ncol(sample_matrix) == 1) {
    aggregated_matrix[, idx] <- sample_matrix[, 1]
  } else {
    aggregated_matrix[, idx] <- rowMeans(sample_matrix, na.rm = TRUE)
  }

  aggregated_matrix[is.nan(aggregated_matrix[, idx]), idx] <- NA_real_

  tenplex_ids <- sample_rows$tenplex_id[!is.na(sample_rows$tenplex_id)]
  tenplex_ids <- sort(unique(tenplex_ids))
  tenplex_label <- if (length(tenplex_ids)) {
    paste(tenplex_ids, collapse = "|")
  } else {
    NA_character_
  }

  tmt_channels <- if ("TMT.Channel" %in% names(sample_rows)) {
    channel_vals <- sample_rows$TMT.Channel
    channel_vals <- channel_vals[!is.na(channel_vals)]
    channel_vals <- sort(unique(channel_vals))
    if (length(channel_vals)) paste(channel_vals, collapse = "|") else NA_character_
  } else {
    NA_character_
  }

  col_metadata_list[[idx]] <- data.table::data.table(
    sampleid = sample_name,
    depmap_id = first_non_missing(sample_rows$depmap_id),
    cell_line = first_non_missing(sample_rows$cell_line),
    tissue = first_non_missing(sample_rows$tissue),
    dataset_sample_id = paste(sample_cols, collapse = "|"),
    n_measurements = nrow(sample_rows),
    plex_id = tenplex_label,
    tmt_channels = tmt_channels
  )
}

col_metadata_dt <- data.table::rbindlist(col_metadata_list, use.names = TRUE, fill = TRUE)

col_metadata_dt[
  ,
  batchid := ifelse(is.na(plex_id) | plex_id == "", NA_character_, plex_id)
]

setcolorder(
  col_metadata_dt,
  c(
    "sampleid",
    "depmap_id",
    "cell_line",
    "tissue",
    "dataset_sample_id",
    "n_measurements",
    "plex_id",
    "batchid",
    "tmt_channels"
  )
)

# SummarizedExperiment ---------------------------------------------------- #

message("Constructing SummarizedExperiment")

colnames(aggregated_matrix) <- col_metadata_dt$sampleid

rowData <- S4Vectors::DataFrame(
  feature_dt,
  row.names = feature_dt$feature_id
)

colData <- S4Vectors::DataFrame(
  sampleid = col_metadata_dt$sampleid,
  depmap_id = col_metadata_dt$depmap_id,
  cell_line = col_metadata_dt$cell_line,
  tissue = col_metadata_dt$tissue,
  batchid = col_metadata_dt$batchid,
  dataset_sample_id = col_metadata_dt$dataset_sample_id,
  n_measurements = col_metadata_dt$n_measurements,
  plex_id = col_metadata_dt$plex_id,
  tmt_channels = col_metadata_dt$tmt_channels,
  row.names = col_metadata_dt$sampleid
)

mass_spec_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = aggregated_matrix),
  rowData = rowData,
  colData = colData
)

SummarizedExperiment::metadata(mass_spec_se) <- list(
  annotation = "proteomics",
  datatype = "proteomics",
  technology = "tmt_mass_spectrometry",
  class = "SummarizedExperiment",
  filename = basename(INPUT$quant),
  release = release,
  data_source = cfg,
  numSamples = ncol(mass_spec_se),
  numFeatures = nrow(mass_spec_se),
  generated_on = as.character(Sys.Date())
)

message("Mass spectrometry SummarizedExperiment summary:")
show(mass_spec_se)

# Outputs ----------------------------------------------------------------- #

dir.create(dirname(OUTPUT$massspec_se), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT$feature_metadata), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT$sample_metadata), recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(proteomics.massspec = mass_spec_se),
  file = OUTPUT$massspec_se
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

message("Saved outputs:")
message(" - SummarizedExperiment list: ", OUTPUT$massspec_se)
message(" - Feature metadata: ", OUTPUT$feature_metadata)
message(" - Sample metadata: ", OUTPUT$sample_metadata)

if (exists("snakemake")) {
  sink()
}
