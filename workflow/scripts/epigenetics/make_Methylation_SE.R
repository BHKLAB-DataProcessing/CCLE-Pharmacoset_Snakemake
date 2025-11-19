## ------------------- Parse Snakemake Object ------------------- ##
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  if (length(snakemake@log)) {
    sink(
      snakemake@log[[1]],
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
  if (file.exists("snapshots/make_Methylation_SE.RData")) {
    load("snapshots/make_Methylation_SE.RData")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(SummarizedExperiment)
  library(GenomicRanges)
})

data.table::setDTthreads(THREADS)

## ------------------- Helper Functions ------------------- ##
normalize_seqname <- function(seqname) {
  seqname <- trimws(seqname)
  seqname <- ifelse(
    grepl("^chr", seqname, ignore.case = FALSE),
    seqname,
    paste0("chr", seqname)
  )
  seqname <- gsub("^chrchr", "chr", seqname, ignore.case = TRUE)
  seqname[seqname == "chrMT"] <- "chrM"
  seqname
}

parse_cpg_sites <- function(cpg_string) {
  coords <- strsplit(cpg_string, ";", fixed = TRUE)
  lapply(coords, function(entries) {
    entries <- trimws(entries[nzchar(entries)])
    chroms <- sub(":.*$", "", entries)
    positions <- suppressWarnings(as.integer(sub("^.*:", "", entries)))
    keep <- !is.na(positions)
    chroms <- chroms[keep]
    positions <- positions[keep]
    if (!length(positions)) {
      return(
        list(
          seqnames = NA_character_,
          start = NA_integer_,
          end = NA_integer_,
          count = 0L
        )
      )
    }
    list(
      seqnames = normalize_seqname(chroms[1]),
      start = min(positions),
      end = max(positions),
      count = length(positions)
    )
  })
}

align_samples <- function(dt, feature_col, valid_samples) {
  sample_cols <- setdiff(
    names(dt),
    c(feature_col, "CpG_sites_hg19", "avg_coverage")
  )
  sample_cols <- sample_cols[!duplicated(sample_cols)]

  missing <- setdiff(sample_cols, valid_samples)
  if (length(missing)) {
    warning(sprintf(
      "Dropping %d samples without metadata for %s. Examples: %s",
      length(missing),
      feature_col,
      paste(utils::head(missing, 10), collapse = ", ")
    ))
  }

  keep_cols <- intersect(sample_cols, valid_samples)
  if (!length(keep_cols)) {
    stop("No overlapping samples between methylation data and metadata.")
  }

  exprs_mat <- as.matrix(dt[, ..keep_cols])
  storage.mode(exprs_mat) <- "double"
  exprs_mat[is.nan(exprs_mat)] <- NA_real_
  rownames(exprs_mat) <- dt[[feature_col]]

  list(
    matrix = exprs_mat,
    samples = keep_cols
  )
}

write_matrix_tsv <- function(mat, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  export_dt <- data.table::as.data.table(mat, keep.rownames = "feature_id")
  data.table::fwrite(
    export_dt,
    file = file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    col.names = TRUE
  )
}

log_and_quit <- function(e) {
  message("ERROR: ", conditionMessage(e))
  tb <- utils::capture.output(sys.calls())
  if (length(tb)) {
    message(paste(tb, collapse = "\n"))
  }
  sink(NULL)
  q(save = "no", status = 1, runLast = FALSE)
}

main <- function() {
  ## ------------------- Load Metadata ------------------- ##
  sample_metadata <- data.table::fread(INPUT$sample_metadata)
  sample_metadata <- sample_metadata[!duplicated(CCLE.sampleid), ]
  valid_samples <- sample_metadata$CCLE.sampleid
  data.table::setkeyv(sample_metadata, "CCLE.sampleid")

  dataset_specs <- list(
    tss_1kb = list(
      file = INPUT$tss_1kb,
      id_col = "locus_id",
      feature_type = "promoter_1kb",
      output = OUTPUT$tss_1kb_matrix
    ),
    tss_cpg_clusters = list(
      file = INPUT$tss_cpg_clusters,
      id_col = "cluster_id",
      feature_type = "promoter_cpg_cluster",
      output = OUTPUT$tss_cpg_clusters_matrix
    ),
    cgi_cpg_clusters = list(
      file = INPUT$cgi_cpg_clusters,
      id_col = "cluster_id",
      feature_type = "cgi_cpg_cluster",
      output = OUTPUT$cgi_cpg_clusters_matrix
    ),
    enhancer_cpg_clusters = list(
      file = INPUT$enhancer_cpg_clusters,
      id_col = "cluster_id",
      feature_type = "enhancer_cpg_cluster",
      output = OUTPUT$enhancer_cpg_clusters_matrix
    )
  )

  dataset_metadata <- snakemake@params$dataset_metadata
  se_outputs <- list()

  for (dataset_name in names(dataset_specs)) {
    spec <- dataset_specs[[dataset_name]]
    message(sprintf("Processing RRBS dataset: %s", dataset_name))

    dt <- data.table::fread(
      spec$file,
      sep = "\t",
      header = TRUE,
      na.strings = c("NA", "NaN", "nan", "")
    )

    dt[[spec$id_col]] <- trimws(dt[[spec$id_col]])
    dt[, CpG_sites_hg19 := trimws(CpG_sites_hg19)]
    dt[, avg_coverage := as.numeric(avg_coverage)]

    if (dataset_name == "tss_1kb") {
      parsed <- stringr::str_match(
        dt[[spec$id_col]],
        "^(.*)_(?:chr)?([0-9A-Za-z]+)_(\\d+)_(\\d+)$"
      )
      gene_symbol <- parsed[, 2]
      chrom <- parsed[, 3]
      start <- suppressWarnings(as.integer(parsed[, 4]))
      end <- suppressWarnings(as.integer(parsed[, 5]))
      seqnames <- normalize_seqname(chrom)
      valid_idx <- !is.na(seqnames) & !is.na(start) & !is.na(end)
      if (any(!valid_idx)) {
        warning(sprintf(
          "Dropping %d promoter loci without valid coordinates.",
          sum(!valid_idx)
        ))
      }
      dt <- dt[valid_idx]
      gene_symbol <- gene_symbol[valid_idx]
      seqnames <- seqnames[valid_idx]
      start <- start[valid_idx]
      end <- end[valid_idx]
      gr <- GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(start = start, end = end),
        strand = "*",
        gene_symbol = gene_symbol,
        locus_id = dt[[spec$id_col]]
      )
      S4Vectors::mcols(gr)$aggregation <- spec$feature_type
      S4Vectors::mcols(gr)$avg_coverage <- dt$avg_coverage
      S4Vectors::mcols(gr)$cpg_sites <- dt$CpG_sites_hg19
    } else if (
      dataset_name %in% c("cgi_cpg_clusters", "enhancer_cpg_clusters")
    ) {
      parsed <- stringr::str_match(
        dt[[spec$id_col]],
        "^(?:chr)?([0-9A-Za-z]+)_(\\d+)_(\\d+)$"
      )
      seqnames <- normalize_seqname(parsed[, 2])
      start <- suppressWarnings(as.integer(parsed[, 3]))
      end <- suppressWarnings(as.integer(parsed[, 4]))
      valid_idx <- !is.na(seqnames) & !is.na(start) & !is.na(end)
      if (any(!valid_idx)) {
        warning(sprintf(
          "Dropping %d CpG clusters without valid coordinates in %s.",
          sum(!valid_idx),
          dataset_name
        ))
      }
      dt <- dt[valid_idx]
      seqnames <- seqnames[valid_idx]
      start <- start[valid_idx]
      end <- end[valid_idx]
      gr <- GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(start = start, end = end),
        strand = "*",
        cluster_id = dt[[spec$id_col]]
      )
      cpg_ranges <- parse_cpg_sites(dt$CpG_sites_hg19)
      S4Vectors::mcols(gr)$num_cpg_sites <- vapply(
        cpg_ranges,
        function(x) x$count,
        integer(1)
      )
      S4Vectors::mcols(gr)$avg_coverage <- dt$avg_coverage
      S4Vectors::mcols(gr)$cpg_sites <- dt$CpG_sites_hg19
    } else if (dataset_name == "tss_cpg_clusters") {
      parsed <- stringr::str_match(
        dt[[spec$id_col]],
        "^(.*)_([0-9]+)$"
      )
      gene_symbol <- parsed[, 2]
      cluster_index <- parsed[, 3]
      cpg_ranges <- parse_cpg_sites(dt$CpG_sites_hg19)
      seqnames <- vapply(cpg_ranges, function(x) x$seqnames, character(1))
      start <- vapply(cpg_ranges, function(x) x$start, integer(1))
      end <- vapply(cpg_ranges, function(x) x$end, integer(1))
      count <- vapply(cpg_ranges, function(x) x$count, integer(1))
      valid_idx <- !is.na(seqnames) & !is.na(start) & !is.na(end)
      if (any(!valid_idx)) {
        warning(sprintf(
          "Dropping %d promoter CpG clusters without coordinates.",
          sum(!valid_idx)
        ))
      }
      dt <- dt[valid_idx]
      seqnames <- seqnames[valid_idx]
      start <- start[valid_idx]
      end <- end[valid_idx]
      gene_symbol <- gene_symbol[valid_idx]
      cluster_index <- cluster_index[valid_idx]
      count <- count[valid_idx]
      gr <- GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(start = start, end = end),
        strand = "*",
        gene_symbol = gene_symbol,
        cluster_label = dt[[spec$id_col]],
        cluster_index = cluster_index
      )
      S4Vectors::mcols(gr)$num_cpg_sites <- count
      S4Vectors::mcols(gr)$avg_coverage <- dt$avg_coverage
      S4Vectors::mcols(gr)$cpg_sites <- dt$CpG_sites_hg19
    } else {
      stop(sprintf("Unhandled dataset: %s", dataset_name))
    }

    sample_alignment <- align_samples(
      dt,
      spec$id_col,
      valid_samples
    )

    coldata <- data.table::data.table(
      sampleid = sample_alignment$samples,
      dataset_sample_id = sample_alignment$samples,
      batchid = rep(NA, length(sample_alignment$samples))
    )

    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        exprs = sample_alignment$matrix
      ),
      rowRanges = gr,
      colData = coldata
    )

    metadata_list <- list(
      annotation = "methylation",
      datatype = dataset_name,
      feature_type = spec$feature_type,
      class = "RangedSummarizedExperiment",
      filename = basename(spec$file),
      data_source = dataset_metadata[[dataset_name]],
      numSamples = ncol(se),
      numFeatures = nrow(se),
      date = Sys.Date()
    )
    se@metadata <- metadata_list

    se_outputs[[paste0("rrbs.", dataset_name)]] <- se
    write_matrix_tsv(sample_alignment$matrix, spec$output)
  }

  message(sprintf(
    "Saving methylation SE list to %s",
    OUTPUT$methylation_se_list
  ))
  dir.create(
    dirname(OUTPUT$methylation_se_list),
    recursive = TRUE,
    showWarnings = FALSE
  )
  saveRDS(se_outputs, OUTPUT$methylation_se_list)
}

tryCatch(
  main(),
  error = log_and_quit
)

sink()
