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
  if (file.exists(file.path("snapshots", "make_Fusion_SE.RData"))) {
    load(file.path("snapshots", "make_Fusion_SE.RData"))
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(S4Vectors)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

fusion_cfg <- PARAMS$dataset_config
preferred_dataset <- fusion_cfg$preferred %||% "filtered"

message("Building fusion SummarizedExperiment")
message(" - Filtered calls: ", INPUT$filtered)
message(" - Unfiltered calls: ", INPUT$unfiltered)

fusion_dt <- data.table::fread(
  if (preferred_dataset == "filtered") INPUT$filtered else INPUT$unfiltered,
  sep = "\t",
  showProgress = FALSE
)

if (!"J_FFPM" %in% names(fusion_dt)) fusion_dt[, J_FFPM := NA_real_]
if (!"S_FFPM" %in% names(fusion_dt)) fusion_dt[, S_FFPM := NA_real_]

required_cols <- c("X.sample", "X.FusionName", "JunctionReadCount", "SpanningFragCount", "BroadID")
missing_cols <- setdiff(required_cols, names(fusion_dt))
if (length(missing_cols)) {
  stop("Fusion table missing required columns: ", paste(missing_cols, collapse = ", "))
}

sample_md <- data.table::fread(INPUT$sample_metadata, showProgress = FALSE)
lookup_cols <- c("CCLE.depMapID", "CCLE.sampleid", "cellosaurus.cellLineName", "CCLE.site_Primary", "CCLE.type")
if (!all(lookup_cols %in% names(sample_md))) {
  stop("Annotated sample metadata missing required columns.", call. = FALSE)
}

sample_lookup <- unique(sample_md[, ..lookup_cols])

fusion_dt[, dataset_sample_id := X.sample]
fusion_dt[, fusion_id := X.FusionName]
fusion_dt[, fusion_uid := make.unique(fusion_id)]

sample_map <- unique(fusion_dt[, .(dataset_sample_id, BroadID)])

setkey(sample_lookup, CCLE.sampleid)
col_meta <- sample_lookup[sample_map, on = c("CCLE.sampleid" = "dataset_sample_id")]

missing_mask <- is.na(col_meta$CCLE.sampleid)
if (any(missing_mask)) {
  setkey(sample_lookup, CCLE.depMapID)
  col_meta[missing_mask, c(
    "CCLE.depMapID",
    "CCLE.sampleid",
    "cellosaurus.cellLineName",
    "CCLE.site_Primary",
    "CCLE.type"
  ) := sample_lookup[
    BroadID,
    .(
      CCLE.depMapID,
      CCLE.sampleid,
      cellosaurus.cellLineName,
      CCLE.site_Primary,
      CCLE.type
    ),
    on = "CCLE.depMapID"
  ]
  ]
}

if (any(is.na(col_meta$CCLE.sampleid))) {
  unresolved <- unique(col_meta$dataset_sample_id[is.na(col_meta$CCLE.sampleid)])
  stop("Failed to map fusion samples to CCLE.sampleid: ",
       paste(unresolved, collapse = ", "))
}

col_meta[, dataset_sample_id := sample_map$dataset_sample_id]
col_meta[, BroadID := sample_map$BroadID]

# Build assay matrices ---------------------------------------------------- #
fusion_ids <- unique(fusion_dt$fusion_uid)
sample_ids <- col_meta$CCLE.sampleid

junction_mat <- matrix(
  0,
  nrow = length(fusion_ids),
  ncol = length(sample_ids),
  dimnames = list(fusion_ids, sample_ids)
)
spanning_mat <- junction_mat
ffpm_j <- junction_mat
ffpm_s <- junction_mat

fusion_dt_subset <- fusion_dt[, .(
  dataset_sample_id,
  fusion_uid,
  JunctionReadCount,
  SpanningFragCount,
  J_FFPM = ifelse("J_FFPM" %in% names(fusion_dt), J_FFPM, NA_real_),
  S_FFPM = ifelse("S_FFPM" %in% names(fusion_dt), S_FFPM, NA_real_)
)]

for (i in seq_len(nrow(fusion_dt_subset))) {
  row <- fusion_dt_subset[i]
  sid <- col_meta$CCLE.sampleid[col_meta$dataset_sample_id == row$dataset_sample_id][1]
  if (is.na(sid)) next
  fid <- row$fusion_uid
  junction_mat[fid, sid] <- row$JunctionReadCount
  spanning_mat[fid, sid] <- row$SpanningFragCount
  ffpm_j[fid, sid] <- row$J_FFPM
  ffpm_s[fid, sid] <- row$S_FFPM
}

# Feature metadata -------------------------------------------------------- #
feature_md <- fusion_dt[, .(
  fusion_id = fusion_uid,
  fusion_name = X.FusionName,
  left_gene = LeftGene,
  right_gene = RightGene,
  left_breakpoint = LeftBreakpoint,
  right_breakpoint = RightBreakpoint,
  splice_type = SpliceType
)] |>
  unique(by = "fusion_id")

rownames_df <- feature_md$fusion_id

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    junction_reads = junction_mat,
    spanning_frags = spanning_mat,
    junction_ffpm = ffpm_j,
    spanning_ffpm = ffpm_s
  ),
  colData = S4Vectors::DataFrame(col_meta[, .(
    sampleid = CCLE.sampleid,
    dataset_sample_id,
    depmap_id = CCLE.depMapID,
    broad_id = BroadID,
    tissue = CCLE.site_Primary,
    type = CCLE.type
  )], row.names = sample_ids),
  rowData = S4Vectors::DataFrame(feature_md, row.names = rownames_df),
  metadata = list(
    annotation = "fusion",
    datatype = preferred_dataset,
    class = "SummarizedExperiment",
    data_source = fusion_cfg,
    date = fusion_cfg[[preferred_dataset]]$release %||% Sys.Date(),
    numSamples = ncol(junction_mat),
    numFeatures = nrow(junction_mat)
  )
)

message("Fusion SE samples: ", ncol(se))
message("Fusion SE features: ", nrow(se))

dir.create(dirname(OUTPUT$fusion_se), showWarnings = FALSE, recursive = TRUE)
saveRDS(se, file = OUTPUT$fusion_se)
data.table::fwrite(feature_md, OUTPUT$feature_metadata, sep = "\t")
data.table::fwrite(as.data.table(colData(se)), OUTPUT$sample_metadata, sep = "\t")

message("Fusion SummarizedExperiment written to ", OUTPUT$fusion_se)

if (length(LOG)) sink()
