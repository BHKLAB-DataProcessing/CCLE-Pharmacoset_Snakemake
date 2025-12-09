## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  LOG <- snakemake@log

  if (length(LOG)) {
    sink(LOG[[1]], split = TRUE, type = c("output", "message"))
  }

  dir.create("snapshots", showWarnings = FALSE, recursive = TRUE)
  save.image(file.path("snapshots", paste0(snakemake@rule, ".RData")))
} else {
  if (file.exists(file.path("snapshots", "make_CNV_segments_RE.RData"))) {
    load(file.path("snapshots", "make_CNV_segments_RE.RData"))
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(RaggedExperiment)
})

message("Reading CNV segments: ", INPUT$seg)
dt <- data.table::fread(INPUT$seg)
if (names(dt)[1] == "") dt[, (1) := NULL]

# Normalize column names from DepMap segment file
setnames(dt, c("ModelID", "CONTIG", "START", "END", "SEGMENT_COPY_NUMBER"),
         c("depmap_id", "chrom", "start", "end", "segment_mean"), skip_absent = TRUE)

required <- c("depmap_id", "chrom", "start", "end", "segment_mean")
missing <- setdiff(required, names(dt))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

dt[, chrom := gsub("^chr", "", chrom)]
gr_list <- split(dt, dt$depmap_id)
gr_list <- lapply(gr_list, function(sub) {
  GRanges(
    seqnames = sub$chrom,
    ranges = IRanges(start = sub$start, end = sub$end),
    segment_mean = sub$segment_mean
  )
})

re <- RaggedExperiment::RaggedExperiment(gr_list)

metadata(re) <- list(annotation = "cnv", datatype = "segments_wgs", source = basename(INPUT$seg))

saveRDS(re, file = OUTPUT$cnv_segments)

if (length(LOG)) sink()
