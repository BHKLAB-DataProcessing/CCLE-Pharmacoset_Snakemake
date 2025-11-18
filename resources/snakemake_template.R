## --------------------------- Setup ---------------------------- ##
rdata_file <-
  gsub("--file=", "", commandArgs()[4]) |>
  basename() |>
  gsub(".R", ".RData", .)

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

  save.image(file.path("resources", rdata_file))
} else {
  load(file.path("resources", rdata_file))
}
