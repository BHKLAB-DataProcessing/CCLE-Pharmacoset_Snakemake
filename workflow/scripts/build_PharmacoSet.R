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

suppressPackageStartupMessages(library(PharmacoGx))
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


# Read in the treatmentResponseExperiment
# ---------------------------------------
message(paste("Loading: ", INPUT$treatmentResponseExperiment, sep = "\n\t"))
tre <- readRDS(INPUT$treatmentResponseExperiment)


# 1.0 Create additional metadata
# ------------------------------
# sampleNames <- c(
#     lapply(MultiAssayExperiment::colnames(mae), unique) |>
#         unlist() |>
#         unique(),
#     colnames(tre) |>
#         unique()
# ) |> unique()
data.table::setkeyv(sampleMetadata, "CCLE.sampleid")
sampleMetadata[, sampleid := cellosaurus.cellLineName]
sampleMetadata <- sampleMetadata[!duplicated(sampleid), ][!is.na(sampleid), ]

sample <- as.data.frame(
  sampleMetadata,
  row.names = sampleMetadata[, sampleid]
)
sample$unique.sampleid <- rownames(sample)

treatmentNames <- c(
  treatmentMetadata[, unique(CCLE.treatmentid)],
  CoreGx::rowData(tre)[, cleanCharacterStrings(unique(treatmentid))]
) |>
  unique()


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
