## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads
  # setup logger if log file is provided
  if (length(snakemake@log) > 0) {
    sink(
      file = snakemake@log[[1]],
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
  if (file.exists("snapshots/build_MultiAssayExperiment.RData")) {
    load("snapshots/build_MultiAssayExperiment.RData")
  }
}
library(MultiAssayExperiment)
library(data.table)

`%||%` <- function(x, y) if (!is.null(x)) x else y


# Read in metadata
# ----------------
message(paste("Loading: ", INPUT$sampleMetadata, sep = "\n\t"))
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)
sampleMetadata <- sampleMetadata[!duplicated(cellosaurus.cellLineName), ]
# Read in the summarized experiments
# ----------------------------------
message(paste("Loading: ", INPUT$summarizedExperiment_lists, sep = "\n\t"))
se_list <- lapply(INPUT$summarizedExperiment_lists, function(x) {
  message(paste("Loading: ", x, sep = "\n\t"))
  readRDS(x)
}) |>
  unlist(recursive = FALSE)

# Ensure each experiment has a unique, non-empty name -------------------- #
name_fallback <- function(se_obj, idx) {
  meta <- tryCatch(se_obj@metadata, error = function(...) list())
  ann <- meta$annotation %||% meta$datatype %||% meta$data_source$dataset %||% NULL
  if (is.null(ann) || ann == "") ann <- paste0("exp", idx)
  ann
}

names(se_list) <- mapply(
  function(nm, obj, idx) {
    if (!is.null(nm) && nm != "") return(nm)
    name_fallback(obj, idx)
  },
  nm = names(se_list),
  obj = se_list,
  idx = seq_along(se_list)
)

# Prefer informative names from metadata$datatype when available
names(se_list) <- mapply(function(obj, nm) {
  dt <- tryCatch(obj@metadata$datatype, error = function(...) NA_character_)
  if (!is.null(dt) && !is.na(dt) && dt != "") return(dt)
  nm
}, obj = se_list, nm = names(se_list))

# TODO: fix this mess...

# Normalize to descriptive, unique assay names
rename_map <- c(
  "genes_tpm" = "rnaseq.gene_tpm",
  "transcripts_tpm" = "rnaseq.transcript_tpm",
  "genes_counts" = "rnaseq.gene_counts",
  "genes_rpkm" = "rnaseq.gene_rpkm",
  "genes" = "cnv.gene_log2",
  "genes.1" = "mutation.gene_binary",
  "proteomics" = "proteomics.rppa",
  "proteomics.1" = "proteomics.massspec_intensity",
  "tss_1kb" = "methylation.tss_1kb",
  "tss_cpg_clusters" = "methylation.tss_cpg_clusters",
  "cgi_cpg_clusters" = "methylation.cgi_cpg_clusters",
  "enhancer_cpg_clusters" = "methylation.enhancer_cpg_clusters",
  "metabolites" = "metabolomics.lcms",
  "mirna_gct" = "mirna.gene_gct",
  "mirna_mimat" = "mirna.mimat_expression",
  "exon_usage_ratio" = "exon.usage_ratio",
  "exon_usage_denominator" = "exon.usage_denominator",
  "segments_wgs" = "cnv.segments_wgs",
  "global_histone_mrm" = "chromatin.histone_mrm",
  "filtered" = "fusion.calls"
)
names(se_list) <- vapply(
  names(se_list),
  function(nm) {
    if (!is.null(rename_map[[nm]])) rename_map[[nm]] else nm
  },
  character(1)
)

if (is.null(names(se_list))) {
  names(se_list) <- paste0("exp", seq_along(se_list))
}

names(se_list)[is.na(names(se_list)) | names(se_list) == ""] <- paste0(
  "exp",
  which(is.na(names(se_list)) | names(se_list) == "")
)

names(se_list) <- make.unique(names(se_list))

# Collect sample IDs across experiments and drop any not in metadata ---- #
data.table::setkeyv(sampleMetadata, "CCLE.sampleid")
sampleMetadata[, sampleid := CCLE.sampleid]
sampleMetadata[, depmap_id := CCLE.depMapID]

map_to_ccle <- function(ids) {
  ids_ccle <- ifelse(
    ids %in% sampleMetadata$sampleid,
    ids,
    sampleMetadata$sampleid[match(ids, sampleMetadata$depmap_id)]
  )
  ids_ccle
}

# Map column names to CCLE sample ids where possible, then drop unmapped
se_list <- lapply(se_list, function(x) {
  mapped <- map_to_ccle(SummarizedExperiment::colnames(x))
  keep <- !is.na(mapped) & mapped != ""
  if (!all(keep)) {
    message("Dropping ", sum(!keep), " unmapped samples for experiment ",
            x@metadata$annotation %||% "")
  }
  if (any(keep)) {
    x <- x[, keep]
    colnames(x) <- mapped[keep]
  } else {
    # preserve empty experiment with zero samples
    colnames(x) <- character(0)
    x <- x[, 0]
  }
  x
})


# Build MultiAssayExperiment
# --------------------------
summarizedExperimentLists <- se_list
summarizedExperimentLists <- lapply(seq_along(summarizedExperimentLists), function(i) {
  x <- summarizedExperimentLists[[i]]
  existing_cd <- as.data.frame(SummarizedExperiment::colData(x))

  if (!nrow(existing_cd)) {
    existing_cd <- data.frame(sampleid = colnames(x))
  }

  original_sampleid <- if ("sampleid" %in% names(existing_cd)) {
    as.character(existing_cd$sampleid)
  } else {
    colnames(x)
  }

  if (!"dataset_sample_id" %in% names(existing_cd)) {
    existing_cd$dataset_sample_id <- original_sampleid
  }

  if (!"ccle_sample_id" %in% names(existing_cd)) {
    existing_cd$ccle_sample_id <- original_sampleid
  }

  # Join expanded sample metadata to add rich columns where available
  meta_join <- sampleMetadata[data.table(sampleid = colnames(x)), on = "sampleid"]
  meta_join[, dataset_sample_id := existing_cd$dataset_sample_id]
  meta_join[, ccle_sample_id := existing_cd$ccle_sample_id]

  # Keep original sample order
  meta_join <- meta_join[match(colnames(x), meta_join$sampleid)]

  # Ensure required columns exist
  if (!"batchid" %in% names(meta_join)) {
    meta_join[, batchid := NA]
  }

  meta_join_df <- as.data.frame(meta_join)
  stopifnot(nrow(meta_join_df) == length(colnames(x)))
  tryCatch({
    rownames(meta_join_df) <- colnames(x)
  }, error = function(e) {
    stop(
      paste0(
        "rownames assignment failed for experiment ",
        names(summarizedExperimentLists)[i], "; nrow=",
        nrow(meta_join_df), " nsamples=", length(colnames(x)),
        " first names: ", paste(head(colnames(x), 3), collapse = ",")
      ),
      call. = FALSE
    )
  })

  x <- tryCatch({
    if (methods::is(x, "RaggedExperiment")) {
      x
    } else {
      x@colData <- MultiAssayExperiment::DataFrame(
        meta_join_df,
        row.names = colnames(x)
      )
      x
    }
  }, error = function(e) {
    stop(
      paste0(
        "colData assignment failed for experiment ",
        names(summarizedExperimentLists)[i],
        " (class=", class(x), "): ncol=", ncol(x),
        " meta_rows=", nrow(meta_join_df), " message=", conditionMessage(e)
      ),
      call. = FALSE
    )
  })
  x
})
names(summarizedExperimentLists) <- names(se_list)
ExpList <- MultiAssayExperiment::ExperimentList(summarizedExperimentLists)
message(paste("ExperimentList:", capture.output(show(ExpList)), sep = "\n\t"))

data.table::setkeyv(sampleMetadata, "sampleid")
# Create a sample map for each experiment in the ExperimentList
sampleMapList <- lapply(summarizedExperimentLists, function(se) {
  data.frame(
    primary = colnames(se),
    colname = colnames(se),
    stringsAsFactors = FALSE
  )
})

names(sampleMapList) <- names(ExpList)
message(paste(
  "Sample map list:",
  capture.output(str(sampleMapList)),
  sep = "\n\t"
))

# Metadata List
# go through each experiment, extract the metadata and add it to a list
metadata_list <- lapply(summarizedExperimentLists, function(se) {
  metadata_ <- slot(se, "metadata")
  if (metadata_$annotation == "rnaseq") {
    metadata_ <- metadata_[-which(names(metadata_) == "sessionInfo")]
  }
  metadata_
})

# create a data frame for coldata including sampleids and batch ids
sampleMetadata <- sampleMetadata[!duplicated(cellosaurus.cellLineName), ]
colData <- as.data.frame(sampleMetadata, row.names = sampleMetadata$sampleid)
rownames(colData) <- sampleMetadata$sampleid

colData$batchid <- 1

message(sprintf(
  "Column data has %d rows and %d columns",
  nrow(colData),
  ncol(colData)
))
str(colData)

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = ExpList,
  colData = colData,
  sampleMap = MultiAssayExperiment::listToMap(sampleMapList),
  metadata = metadata_list
)


# Write Output
# ------------
message("Saving MultiAssayExperiment to: ", OUTPUT$mae)
saveRDS(mae, file = OUTPUT$mae)
sink()
