options(
  repos = c(CRAN = Sys.getenv("CRAN_REPO", "https://cloud.r-project.org"))
)

current_r_minor <- paste(
  R.version$major,
  strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1],
  sep = "."
)
pak_desc <- tryCatch(
  utils::packageDescription("pak"),
  error = function(e) NA
)
pak_built <- if (is.na(pak_desc[1])) NA_character_ else pak_desc[["Built"]]
pak_built_minor <- sub("^R ([0-9]+\\.[0-9]+).*$", "\\1", pak_built)
pak_needs_install <- is.na(pak_built_minor) ||
  !identical(pak_built_minor, current_r_minor)

if (pak_needs_install) {
  message("Installing pak from CRAN")
  install.packages("pak")
}

github_spec <- function(repo, envvar, default_ref = "") {
  ref <- Sys.getenv(envvar, unset = default_ref)

  if (nzchar(ref)) {
    paste0(repo, "@", ref)
  } else {
    repo
  }
}

cran_packages <- c(
  "BiocManager",
  "DT",
  "Polychrome",
  "UpSetR",
  "data.table",
  "dplyr",
  "forcats",
  "ggplot2",
  "ggridges",
  "ggsci",
  "here",
  "httr",
  "jsonlite",
  "knitr",
  "matrixStats",
  "paletteer",
  "patchwork",
  "pheatmap",
  "purrr",
  "readr",
  "readxl",
  "rlang",
  "rmarkdown",
  "scales",
  "scico",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse"
)

bioc_packages <- c(
  "AnnotationDbi",
  "BiocGenerics",
  "BiocParallel",
  "GenomicRanges",
  "IRanges",
  "MultiAssayExperiment",
  "RaggedExperiment",
  "S4Vectors",
  "SummarizedExperiment",
  "org.Hs.eg.db",
  "rtracklayer"
)

github_packages <- c(
  github_spec("bhklab/CoreGx", "COREGX_REF"),
  github_spec("bhklab/PharmacoGx", "PHARMACOGX_REF"),
  github_spec("bhklab/AnnotationGx", "ANNOTATIONGX_REF", "v0.0.0.9097")
)

requested_packages <- c(
  cran_packages,
  paste0("bioc::", bioc_packages),
  github_packages
)

message("Installing ", length(requested_packages), " R package specs with pak")
pak::pkg_install(
  requested_packages,
  upgrade = FALSE,
  ask = FALSE,
  dependencies = c("Depends", "Imports", "LinkingTo")
)
