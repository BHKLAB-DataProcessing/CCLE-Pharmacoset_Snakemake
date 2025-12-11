# Shared palette helpers for all QC notebooks

suppressPackageStartupMessages({
  if (!requireNamespace("paletteer", quietly = TRUE)) {
    stop("Package 'paletteer' is required for QC palettes; please install it.")
  }
  library(paletteer)
})

# Base categorical palette (20-color D3 Category20) used throughout QC
palette_classic <- as.character(paletteer::paletteer_d("ggsci::category20_d3"))

# Extended palette for tissues; remove grays to keep categories vibrant
tissue_palette_name <- "Polychrome::green_armytage"
palette_extended_all <- as.character(paletteer::paletteer_d(tissue_palette_name))
strip_alpha <- function(hex) substr(hex, 1, 7)
is_gray <- function(hex) {
  rgb <- grDevices::col2rgb(strip_alpha(hex))
  rgb[1, ] == rgb[2, ] & rgb[2, ] == rgb[3, ]
}
palette_extended <- palette_extended_all[!is_gray(palette_extended_all)]
na_color <- "#B0B0B0FF"

# Simple helper to repeat the base palette to requested length
fill_palette <- function(n) rep(palette_classic, length.out = n)

# Lazily build a tissue palette that is consistent across notebooks
tissue_color_map <- NULL
tissue_palette <- function(tissues) {
  available <- tissues[!is.na(tissues)]
  available <- as.character(available)
  available <- ifelse(available == "", "unspecified", available)
  if (exists("tissue_color_map", inherits = FALSE) && length(tissue_color_map)) {
    cols <- tissue_color_map
  } else {
    unique_tissues <- sort(unique(available))
    cols <- setNames(
      rep(palette_extended, length.out = length(unique_tissues)),
      unique_tissues
    )
  }
  present <- cols[names(cols) %in% unique(available)]
  if (any(is.na(tissues))) {
    present <- c(present, `NA` = na_color)
  }
  present
}

# Convenient short-hand colors to keep accents aligned across reports
col1 <- palette_classic[1]
col2 <- palette_classic[2]
col3 <- palette_classic[3]
col4 <- palette_classic[4]
col5 <- palette_classic[5]
col6 <- palette_classic[6]
col7 <- palette_classic[7]
col8 <- palette_classic[8]
