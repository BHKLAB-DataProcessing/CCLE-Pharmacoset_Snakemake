# Packages used in this script:
# - data.table
# - snakemake
# - rtracklayer
# - GenomicRanges
# - IRanges
# - S4Vectors

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
}


# oncomapAssay <- data.table::fread(INPUT$oncomapAssay)

# -- 1. Hybrid Capture

# Code taken from https://github.com/BHKLAB-DataProcessing/get_CCLE/blob/master/get_CCLEP.R
#  mostly unchanged, except for commented out lines

hybridCapture <- read.csv(INPUT$hybridCapture, sep = "\t")
# mut <- hybridCapture[ , c("Hugo_Symbol", "Tumor_Sample_Barcode", "Protein_Change"), drop=FALSE]
mut <- hybridCapture
mut[!is.na(mut) & mut == ""] <- NA
mut[
  is.na(mut[, "Protein_Change"]) | mut[, "Protein_Change"] == "",
  "Protein_Change"
] <- "wt"
mut[
  !is.na(mut[, "Protein_Change"]) &
    (mut[, "Protein_Change"] == "p.?" | mut[, "Protein_Change"] == "p.0?"),
  "Protein_Change"
] <- NA
# mut <- mut[complete.cases(mut), , drop=FALSE]
myx <- !duplicated(paste(
  mut[, c("Tumor_Sample_Barcode")],
  mut[, c("Hugo_Symbol")],
  mut[, c("Protein_Change")],
  sep = "///"
))
mut <- mut[myx, , drop = FALSE]


# -- 2. Oncomap
# Code taken from https://github.com/BHKLAB-DataProcessing/get_CCLE/blob/master/get_CCLEP.R
#  mostly unchanged, except for commented out lines
oncomap <- read.csv(INPUT$oncomap, sep = "\t")
# mut2 <- oncomap[ , c("Hugo_Symbol", "Tumor_Sample_Barcode", "Protein_Change"), drop=FALSE]
mut2 <- oncomap

mut2[!is.na(mut2) & mut2 == ""] <- NA
# mut2[is.na(mut2[ , "Protein_Change"]) | mut2[ , "Protein_Change"] == "", "Protein_Change"] <- "wt"
mut2[
  !is.na(mut2[, "Protein_Change"]) &
    (mut2[, "Protein_Change"] == "p.?" | mut2[, "Protein_Change"] == "p.0?"),
  "Protein_Change"
] <- NA
# mut2 <- mut2[complete.cases(mut2), , drop=FALSE]
myx <- !duplicated(paste(
  mut2[, c("Tumor_Sample_Barcode")],
  mut2[, c("Hugo_Symbol")],
  mut2[, c("Protein_Change")],
  sep = "///"
))
mut2 <- mut2[myx, , drop = FALSE]


# -- 3. Combine
# mutation <- rbind(mut, mut2)
# combine mut and mut2 and only keep columns that are in both
commonCols <- intersect(colnames(mut), colnames(mut2))

mutation <- rbind(mut[, commonCols], mut2[, commonCols])

ucell <- sort(unique(mutation[, "Tumor_Sample_Barcode"]))
ugene <- sort(unique(mutation[, "Hugo_Symbol"]))
dd <- matrix(
  "wt",
  nrow = length(ucell),
  ncol = length(ugene),
  dimnames = list(ucell, ugene)
)

mm <- 1:nrow(mutation)
ff <- TRUE
while (length(mm) > 1) {
  myx <- !duplicated(paste(
    mutation[mm, c("Tumor_Sample_Barcode")],
    mutation[mm, c("Hugo_Symbol")],
    sep = "///"
  ))
  if (ff) {
    dd[as.matrix(mutation[
      mm[myx],
      c("Tumor_Sample_Barcode", "Hugo_Symbol")
    ])] <- mutation[mm[myx], "Protein_Change"]
    ff <- FALSE
  } else {
    dd[as.matrix(mutation[
      mm[myx],
      c("Tumor_Sample_Barcode", "Hugo_Symbol")
    ])] <- paste(
      dd[as.matrix(mutation[
        mm[myx],
        c("Tumor_Sample_Barcode", "Hugo_Symbol")
      ])],
      mutation[mm[myx], "Protein_Change"],
      sep = "///"
    )
  }
  mm <- mm[!myx]
}
## check for inconsistencies (wt + mutations)
iix <- grep("///", dd)
for (iii in iix) {
  x <- sort(unique(unlist(strsplit(dd[iii], split = "///"))))
  if (length(x) > 1) {
    x <- x[!is.element(x, "wt")]
  }
  dd[iii] <- paste(x, collapse = "///")
}


# nn <- sampleinfo[match(rownames(dd), sampleinfo[ , "CCLE.name"]), "cellid"]

# ## remove if we do not have cell line identifier
# dd <- dd[!is.na(nn), , drop=FALSE]
# #nn <- as.character(matchToIDTable(ids=nn, tbl=curationCell, column = "CCLE.cellid", returnColumn = "unique.cellid"))
# rownames(dd) <- nn[!is.na(nn)]
mut_assay <- t(dd)


out <- list(
  assay = mut_assay,
  data = mutation
)

saveRDS(out, file = OUTPUT$preprocessedMutation)
