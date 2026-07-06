args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: query_chembl_mechanism.R <chembl_id> <output_rds>",
    call. = FALSE
  )
}

chembl_id <- args[[1]]
output_rds <- args[[2]]

result <- AnnotationGx::getChemblMechanism(chembl_id)
result <- data.table::as.data.table(result)
saveRDS(result, output_rds)
