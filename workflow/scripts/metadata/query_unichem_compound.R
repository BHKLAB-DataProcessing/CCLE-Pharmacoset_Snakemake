args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: query_unichem_compound.R <compound> <source_id> <output_rds>",
    call. = FALSE
  )
}

compound <- args[[1]]
source_id <- as.integer(args[[2]])
output_rds <- args[[3]]

result <- AnnotationGx::queryUnichemCompound(
  type = "sourceID",
  compound = compound,
  sourceID = source_id
)
saveRDS(result, output_rds)
