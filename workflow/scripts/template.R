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
} else {
  if (file.exists(file.path("snapshots/", paste0(snakemake@rule, ".RData")))) {
    load(file.path("snapshots/", paste0(snakemake@rule, ".RData")))
  }
}
