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
}
if (file.exists("snapshots/build_PharmacoSet.RData")) {
  load("snapshots/build_PharmacoSet.RData")
}

suppressPackageStartupMessages({
  library(PharmacoGx)
  library(MultiAssayExperiment)
})
snakemake@source("metadata/cleanCharacterStrings.R") # for cleanCharacterStrings()

# Read in metadata
# ----------------
message(paste("Loading: ", INPUT$sampleMetadata, sep = "\n\t"))
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)

message(paste("Loading: ", INPUT$treatmentMetadata, sep = "\n\t"))
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)

# Read in the MAE
# ---------------
message(paste("Loading: ", INPUT$multiAssayExperiment, sep = "\n\t"))
mae <- readRDS(INPUT$multiAssayExperiment)

# PharmacoGx PharmacoSet2 expects SummarizedExperiment assays only
keep_exps <- vapply(
  experiments(mae),
  function(x) methods::is(x, "SummarizedExperiment"),
  logical(1)
)
if (!all(keep_exps)) {
  message("Dropping non-SummarizedExperiment assays before PharmacoSet build: ",
          paste(names(experiments(mae))[!keep_exps], collapse = ", "))
}
keep_names <- names(experiments(mae))[keep_exps]
exp_list <- experiments(mae)[keep_names]
sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
sm <- sm[sm$assay %in% keep_names, , drop = FALSE]
sm$assay <- as.character(sm$assay)

# Standardize descriptive assay names to match MAE build
# We support both the legacy names (rna/rrbs/etc) and the current MAE names.
rename_map <- c(
  # Legacy -> new
  "rnaseq.transcripts_tpm" = "rnaseq.transcript_tpm",
  "rnaseq.genes_tpm" = "rnaseq.gene_tpm",
  "rnaseq.genes_counts" = "rnaseq.gene_counts",
  "rnaseq.genes_rpkm" = "rnaseq.gene_rpkm",
  "cnv.genes" = "cnv.gene_log2",
  "mut.genes" = "mutation.gene_binary",
  "proteomics.rppa" = "proteomics.rppa",
  "proteomics.massspec" = "proteomics.massspec_intensity",
  "rrbs.tss_1kb" = "methylation.tss_1kb",
  "rrbs.tss_cpg_clusters" = "methylation.tss_cpg_clusters",
  "rrbs.cgi_cpg_clusters" = "methylation.cgi_cpg_clusters",
  "rrbs.enhancer_cpg_clusters" = "methylation.enhancer_cpg_clusters",
  "metabolomics.lcms" = "metabolomics.lcms",
  "mirna" = "mirna.gene_gct",
  "mirna.1" = "mirna.mimat_expression",
  "exon" = "exon.usage_ratio",
  "exon.1" = "exon.usage_denominator",
  "chromatin" = "chromatin.histone_mrm",
  "fusion" = "fusion.calls",
  # Already-standard names -> themselves (idempotent)
  "rnaseq.transcript_tpm" = "rnaseq.transcript_tpm",
  "rnaseq.gene_tpm" = "rnaseq.gene_tpm",
  "rnaseq.gene_counts" = "rnaseq.gene_counts",
  "rnaseq.gene_rpkm" = "rnaseq.gene_rpkm",
  "cnv.gene_log2" = "cnv.gene_log2",
  "cnv.gene_log2.1" = "mutation.gene_binary",
  "mutation.gene_binary" = "mutation.gene_binary",
  "proteomics.massspec_intensity" = "proteomics.massspec_intensity",
  "proteomics.rppa.1" = "proteomics.massspec_intensity",
  "methylation.tss_1kb" = "methylation.tss_1kb",
  "methylation.tss_cpg_clusters" = "methylation.tss_cpg_clusters",
  "methylation.cgi_cpg_clusters" = "methylation.cgi_cpg_clusters",
  "methylation.enhancer_cpg_clusters" = "methylation.enhancer_cpg_clusters",
  "mirna.gene_gct" = "mirna.gene_gct",
  "mirna.mimat_expression" = "mirna.mimat_expression",
  "exon.usage_ratio" = "exon.usage_ratio",
  "exon.usage_denominator" = "exon.usage_denominator",
  "cnv.segments_wgs" = "cnv.segments_wgs",
  "chromatin.histone_mrm" = "chromatin.histone_mrm",
  "fusion.calls" = "fusion.calls"
)

# Apply rename to experiments and sampleMap without triggering harmonize drops
new_exp_names <- make.unique(vapply(
  names(exp_list),
  function(n) if (!is.null(rename_map[[n]])) rename_map[[n]] else n,
  character(1)
))
names(exp_list) <- new_exp_names

sm$assay <- vapply(
  sm$assay,
  function(a) if (!is.null(rename_map[[a]])) rename_map[[a]] else a,
  character(1)
)

mae <- MultiAssayExperiment(
  experiments = exp_list,
  colData = colData(mae),
  sampleMap = sm,
  metadata = metadata(mae)
)

# Read in the treatmentResponseExperiment
# ---------------------------------------
message(paste("Loading: ", INPUT$treatmentResponseExperiment, sep = "\n\t"))
tre <- readRDS(INPUT$treatmentResponseExperiment)


# 1.0 Create additional metadata
# ------------------------------
data.table::setkeyv(sampleMetadata, "CCLE.sampleid")
sampleMetadata[, sampleid := CCLE.sampleid]
sampleMetadata <- sampleMetadata[!duplicated(sampleid), ][!is.na(sampleid), ]

sample <- as.data.frame(
  sampleMetadata,
  row.names = sampleMetadata[, sampleid]
)
sample$unique.sampleid <- rownames(sample)

data.table::setkeyv(treatmentMetadata, "CCLE.treatmentid")
treatmentMetadata <- treatmentMetadata[!duplicated(CCLE.treatmentid), ]
treatmentMetadata[, treatmentid := CCLE.treatmentid]

treatment <- data.frame(
  treatmentMetadata,
  row.names = treatmentMetadata[, CCLE.treatmentid]
)
treatment$unique.treatmentid <- rownames(treatment)


# _.0 Build PharmacoSet
# ---------------------

name <- "CCLE_(2019)"

pset <- PharmacoGx::PharmacoSet2(
  name = name,
  treatment = treatment,
  sample = sample,
  molecularProfiles = mae,
  treatmentResponse = tre,
  perturbation = list(),
  curation = list(
    sample = sample,
    treatment = treatment,
    tissue = data.frame()
  ),
  datasetType = "sensitivity"
)

# Preserve readable sample IDs in `unique.sampleid`
pset@sample$unique.sampleid <- pset@sample$sampleid
pset@curation$sample$unique.sampleid <- pset@curation$sample$sampleid

message(paste(capture.output(show(pset)), collapse = "\n\t"))

message("Object Size (pset):")
object.size(pset) |> print(unit = "auto")
print(paste("Saving PharmacoSet object to", OUTPUT[[1]]))
dir.create(dirname(OUTPUT[[1]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(pset, file = OUTPUT[[1]])
