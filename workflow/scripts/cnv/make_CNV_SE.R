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

  save.image(
    file.path("resources/", paste0(snakemake@rule, ".RData"))
  )
} else {
  load(file.path("resources/make_CNV_SE.RData"))
}
## ------------------- LOAD INPUT ------------------- ##

message("Loading CNV data")
cnvDt <- data.table::fread(
  INPUT$cnv,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

cnv_seqnames <- paste0("chr", cnvDt$CHR)

message("Creating GRanges object for CNV data")
cnv_gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(cnv_seqnames),
  ranges = IRanges::IRanges(start = cnvDt$CHRLOC, end = cnvDt$CHRLOCEND),
  strand = S4Vectors::Rle("*") # Assuming strand information is not available in your CNV data
)
print(cnv_gr)

message("Creating SummarizedExperiment object for CNV data")
cnv_matrix <- as.matrix(
  cnvDt[, -c("EGID", "SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"), with = FALSE],
  rownames = cnvDt[["SYMBOL"]]
)

cnv_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    exprs = cnv_matrix
  ),
  rowRanges = cnv_gr,
  colData = data.table::data.table(
    sampleid = colnames(cnv_matrix),
    # make a column called batchid that is full of NAs
    batchid = rep(NA, ncol(cnv_matrix))
  )
)

message("Adding metadata to CNV SummarizedExperiment object")
(cnv_se@metadata <- list(
  annotation = "cnv",
  datatype = "genes",
  class = "RangedSummarizedExperiment",
  filename = "CCLE_copynumber_byGene_2013-12-03.txt",
  data_source = snakemake@config$molecularProfiles$cnv,
  date = Sys.Date(),
  numSamples = ncol(cnv_se),
  numGenes = nrow(cnv_se)
))

print(cnv_se)


## ------------------- SAVE OUTPUT ------------------- ##
message("Saving CNV SummarizedExperiment object")
se_list <- list(
  cnv.genes = cnv_se
)

saveRDS(se_list, file = OUTPUT$CNV_SE)
