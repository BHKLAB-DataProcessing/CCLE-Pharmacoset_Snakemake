## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  if (length(snakemake@log) > 0) {
    sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
  }

  dir.create("snapshots", showWarnings = FALSE, recursive = TRUE)
  save.image(
    file.path("snapshots", paste0(snakemake@rule, ".RData"))
  )
} else {
  if (file.exists(file.path("snapshots", "make_RPPA_SE.RData"))) {
    load(file.path("snapshots", "make_RPPA_SE.RData"))
  }
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
dataset_name <- "tcpa_rppa500"
rppa <- rppa_cfg[[dataset_name]]

if (is.null(rppa)) {
  stop("Configured RPPA dataset 'tcpa_rppa500' is missing.", call. = FALSE)
}

release <- rppa_cfg$release %||% rppa$release %||% NA_character_

data_source <- rppa
data_source$release <- data_source$release %||% release
data_source$dataset <- dataset_name

message("Starting RPPA SummarizedExperiment construction")
message("Using RPPA dataset: ", dataset_name)

## ------------------- Load Inputs ------------------- ##
message("Reading RPPA matrix: ", INPUT$matrix)
rppa_dt <- data.table::fread(INPUT$matrix, showProgress = FALSE)

if (!"sample_id" %in% names(rppa_dt)) {
  stop("RPPA matrix is missing 'sample_id' column.", call. = FALSE)
}

message("Rows (samples): ", nrow(rppa_dt))
message("Columns (features + sample id): ", ncol(rppa_dt))

message("Reading annotated sample metadata: ", INPUT$sample_metadata)
sample_metadata <- data.table::fread(
  INPUT$sample_metadata,
  showProgress = FALSE
)

if (
  !all(
    c("CCLE.depMapID", "CCLE.sampleid", "cellosaurus.cellLineName") %in%
      colnames(sample_metadata)
  )
) {
  stop("Annotated sample metadata is missing required columns.", call. = FALSE)
}

sample_lookup <- unique(
  sample_metadata[, .(
    CCLE.depMapID,
    CCLE.sampleid,
    cellosaurus.cellLineName,
    CCLE.site_Primary,
    CCLE.type
  )],
  by = "CCLE.depMapID"
)

sample_lookup_sample <- data.table::copy(sample_lookup)
setkey(sample_lookup_sample, CCLE.sampleid)

dataset_sample_ids <- rppa_dt$sample_id

col_meta <- sample_lookup_sample[
  data.table::data.table(CCLE.sampleid = dataset_sample_ids),
  on = "CCLE.sampleid"
]

col_meta[, dataset_sample_id := dataset_sample_ids]

missing_mask <- is.na(col_meta$CCLE.depMapID)

if (any(missing_mask)) {
  replacements <- data.table::data.table(
    dataset_sample_id = c("KE97_STOMACH", "NCIH684_LARGE_INTESTINE"),
    CCLE.sampleid = c(
      "KE97_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
      "NCIH684_LIVER"
    )
  )

  replacement_hits <- replacements[
    data.table::data.table(
      dataset_sample_id = dataset_sample_ids[missing_mask]
    ),
    on = "dataset_sample_id",
    nomatch = 0L
  ]

  if (nrow(replacement_hits)) {
    replacement_lookup <- sample_lookup_sample[
      replacement_hits$CCLE.sampleid,
      nomatch = 0L
    ]

    col_meta[
      dataset_sample_id %in% replacement_hits$dataset_sample_id,
      c(
        "CCLE.depMapID",
        "CCLE.sampleid",
        "cellosaurus.cellLineName",
        "CCLE.site_Primary",
        "CCLE.type"
      ) :=
        replacement_lookup[
          match(
            dataset_sample_id,
            replacement_hits$dataset_sample_id
          ),
          .(
            CCLE.depMapID,
            CCLE.sampleid,
            cellosaurus.cellLineName,
            CCLE.site_Primary,
            CCLE.type
          )
        ]
    ]

    message("Applied hardcoded mapping for known TCPA sample identifiers.")
  }

  missing_mask <- is.na(col_meta$CCLE.depMapID)

  if (any(missing_mask)) {
    unresolved <- unique(col_meta$dataset_sample_id[missing_mask])
    stop(
      "Failed to map RPPA samples to CCLE identifiers: ",
      paste(unresolved, sep = ", ", collapse = ", ")
    )
  }
}

depmap_ids <- col_meta$CCLE.depMapID
sample_ids <- col_meta$CCLE.sampleid

# Build assay matrix (features x samples)
rppa_value_cols <- setdiff(names(rppa_dt), "sample_id")

if (!length(rppa_value_cols)) {
  stop(
    "RPPA matrix does not contain feature columns after removing identifiers.",
    call. = FALSE
  )
}

rppa_values <- as.matrix(rppa_dt[, ..rppa_value_cols])
rownames(rppa_values) <- dataset_sample_ids
rppa_matrix <- t(rppa_values)

if (!is.numeric(rppa_matrix[1, 1])) {
  stop("RPPA matrix failed to convert to numeric values.", call. = FALSE)
}

colnames(rppa_matrix) <- sample_ids
col_meta[, dataset_sample_id := dataset_sample_ids]

extract_tissue_slug <- function(ids) {
  out <- ifelse(
    !is.na(ids) & grepl("_", ids),
    sub("^[^_]+_", "", ids),
    NA_character_
  )
  tolower(out)
}

col_meta[, tissue_label := CCLE.site_Primary]

col_meta[
  is.na(tissue_label) | tissue_label == "",
  tissue_label := extract_tissue_slug(CCLE.sampleid)
]

col_meta[
  (is.na(tissue_label) | tissue_label == "") & !is.na(dataset_sample_id),
  tissue_label := extract_tissue_slug(dataset_sample_id)
]

col_meta[
  (is.na(tissue_label) | tissue_label == "") & !is.na(CCLE.type) & CCLE.type != "",
  tissue_label := CCLE.type
]

col_meta[
  is.na(tissue_label) | tissue_label == "",
  tissue_label := "unknown"
]

# Row / feature metadata -------------------------------------------------- #
feature_ids <- rownames(rppa_matrix)

feature_metadata <- data.table::data.table(
  feature_id = feature_ids
)

feature_metadata[,
  clean_id := gsub("^X", "", feature_id)
]

feature_metadata[,
  base_id := gsub("(P[STY][0-9].*)$", "", clean_id)
]

feature_metadata[
  base_id == "",
  base_id := clean_id
]

feature_metadata[,
  is_phospho := grepl("P[STY][0-9]", clean_id)
]

feature_metadata[,
  phosphorylation_site := ifelse(
    is_phospho,
    sub("^(.*?)(P[STY][0-9].*)$", "\\2", clean_id),
    NA_character_
  )
]

feature_metadata[,
  target_hint := sub("(P[STY][0-9].*)$", "", clean_id)
]

feature_metadata[
  target_hint == "",
  target_hint := clean_id
]

alias_keys <- unique(c(
  feature_metadata$clean_id,
  feature_metadata$base_id,
  feature_metadata$target_hint
))

alias_map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(alias_keys),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "ALIAS"
)

alias_map_dt <- if (is.null(alias_map)) {
  data.table::data.table(ALIAS = character())
} else {
  data.table::as.data.table(alias_map)
}

annotations_clean <- alias_map_dt[
  feature_metadata,
  on = c("ALIAS" = "clean_id"),
  allow.cartesian = TRUE
][,
  mapping_strategy := "clean"
]

annotations_base <- alias_map_dt[
  feature_metadata,
  on = c("ALIAS" = "base_id"),
  allow.cartesian = TRUE
][,
  mapping_strategy := "base"
]

annotations_hint <- alias_map_dt[
  feature_metadata,
  on = c("ALIAS" = "target_hint"),
  allow.cartesian = TRUE
][,
  mapping_strategy := "target_hint"
]

feature_annotation_matches <- data.table::rbindlist(
  list(annotations_clean, annotations_base, annotations_hint),
  fill = TRUE,
  use.names = TRUE
)

feature_summary <- feature_annotation_matches[
  !is.na(SYMBOL),
  .(
    gene_symbol = paste(sort(unique(na.omit(SYMBOL))), collapse = "|"),
    gene_name = paste(sort(unique(na.omit(GENENAME))), collapse = "|"),
    entrez_id = paste(sort(unique(na.omit(ENTREZID))), collapse = "|"),
    alias_source = paste(sort(unique(na.omit(ALIAS))), collapse = "|"),
    n_mapped_genes = uniqueN(na.omit(SYMBOL)),
    mapping_strategies = paste(sort(unique(mapping_strategy)), collapse = "|")
  ),
  by = feature_id
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
  feature_summary[
    alias_source == "",
    alias_source := NA_character_
  ]
  feature_summary[
    mapping_strategies == "",
    mapping_strategies := NA_character_
  ]
}

feature_metadata <- merge(
  feature_metadata,
  feature_summary,
  by = "feature_id",
  all.x = TRUE
)

feature_metadata[
  is.na(n_mapped_genes),
  n_mapped_genes := 0L
]

setorder_lookup <- match(feature_ids, feature_metadata$feature_id)

if (anyNA(setorder_lookup)) {
  stop(
    "Feature metadata rows could not be aligned with assay feature identifiers.",
    call. = FALSE
  )
}

feature_metadata <- feature_metadata[setorder_lookup]

setcolorder(
  feature_metadata,
  c(
    "feature_id",
    "clean_id",
    "base_id",
    "target_hint",
    "is_phospho",
    "phosphorylation_site",
    "gene_symbol",
    "gene_name",
    "entrez_id",
    "alias_source",
    "mapping_strategies",
    "n_mapped_genes"
  )
)

# Column metadata --------------------------------------------------------- #

col_metadata <- data.table(
  sampleid = sample_ids,
  depmap_id = col_meta$CCLE.depMapID,
  cell_line = col_meta$cellosaurus.cellLineName,
  tissue = col_meta$tissue_label,
  dataset_sample_id = col_meta$dataset_sample_id %||% sample_ids
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
  rowData = S4Vectors::DataFrame(
    feature_metadata,
    row.names = feature_metadata$feature_id
  ),
  colData = S4Vectors::DataFrame(
    sampleid = col_metadata$sampleid,
    depmap_id = col_metadata$depmap_id,
    cell_line = col_metadata$cell_line,
    tissue = col_metadata$tissue,
    batchid = rep(NA, nrow(col_metadata)),
    dataset_sample_id = col_metadata$dataset_sample_id,
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
dir.create(
  dirname(OUTPUT$feature_metadata),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(OUTPUT$sample_metadata),
  recursive = TRUE,
  showWarnings = FALSE
)

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
