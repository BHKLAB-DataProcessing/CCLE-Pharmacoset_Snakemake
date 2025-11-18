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

library(data.table)


message("Loading somatic mutation data")
somatic <- fread(
  INPUT$somatic,
  header = TRUE,
  stringsAsFactors = FALSE
)

if ("Unnamed: 0" %in% names(somatic)) {
  somatic[, "Unnamed: 0" := NULL]
}

metadata_cols <- intersect(
  c(
    "SequencingID",
    "ModelID",
    "ModelConditionID",
    "IsDefaultEntryForModel",
    "IsDefaultEntryForMC"
  ),
  names(somatic)
)

message("Filtering to default model entries")
somatic[, IsDefaultEntryForModel := tolower(IsDefaultEntryForModel)]
somatic[, IsDefaultEntryForMC := tolower(IsDefaultEntryForMC)]
somatic_keep <- somatic[
  IsDefaultEntryForModel == "yes" &
    (is.na(IsDefaultEntryForMC) | IsDefaultEntryForMC == "yes")
]
if (!nrow(somatic_keep)) {
  warning(
    "No rows marked as default for both model and condition; using all rows"
  )
  somatic_keep <- somatic
}

message("Loading sample metadata for ModelID to CCLE.sampleid mapping")
sampleMetadata <- fread(INPUT$sampleMetadata)
sample_lookup <- unique(sampleMetadata[!is.na(CCLE.depMapID), .(
  ModelID = CCLE.depMapID,
  CCLE.sampleid
)])
setkey(sample_lookup, ModelID)

somatic_keep[, sampleid := sample_lookup[ModelID, CCLE.sampleid]]
missing_samples <- somatic_keep[is.na(sampleid), unique(ModelID)]
if (length(missing_samples)) {
  preview <- head(missing_samples, 20)
  warning(sprintf(
    "Using ModelID as sample name for %d samples with no CCLE.sampleid match (showing up to 20): %s",
    length(missing_samples),
    paste(preview, collapse = ", ")
  ))
  somatic_keep[is.na(sampleid), sampleid := ModelID]
}

message("Preparing per-gene mutation calls")
somatic_keep[, HugoSymbol := trimws(HugoSymbol)]
somatic_keep[, ProteinChange := trimws(ProteinChange)]
somatic_keep[is.na(ProteinChange) | ProteinChange == "", ProteinChange := "wt"]

somatic_keep[, gene_call := ifelse(
  ProteinChange == "wt",
  "wt",
  ProteinChange
)]

somatic_keep[
  ,
  gene_call := paste0(
    gene_call,
    ifelse(
      is.na(VepImpact),
      "",
      paste0("(", VepImpact, ")")
    )
  )
]

gene_levels <- sort(unique(somatic_keep$HugoSymbol))
sample_levels <- sort(unique(somatic_keep$sampleid))

mat <- matrix(
  "wt",
  nrow = length(gene_levels),
  ncol = length(sample_levels),
  dimnames = list(gene_levels, sample_levels)
)

somatic_keep <- somatic_keep[
  !is.na(HugoSymbol) & HugoSymbol != "" &
    !is.na(sampleid) & sampleid != ""
]

setorder(somatic_keep, sampleid, HugoSymbol)
somatic_keep[, gene_call := trimws(gene_call)]

collapsed <- somatic_keep[
  ,
  .(gene_call = {
    vals <- unique(gene_call)
    vals <- vals[vals != "wt"]
    if (!length(vals)) "wt" else paste(vals, collapse = "///")
  }),
  by = .(HugoSymbol, sampleid)
]

wide <- data.table::dcast(
  collapsed,
  HugoSymbol ~ sampleid,
  value.var = "gene_call",
  fill = "wt"
)

message("Assembling row (gene) metadata")
gene_meta <- unique(somatic_keep[, .(
  HugoSymbol,
  EntrezGeneID,
  EnsemblGeneID
)])[order(HugoSymbol)]
data.table::setkey(gene_meta, HugoSymbol)

gene_meta <- gene_meta[match(wide$HugoSymbol, HugoSymbol)]

assay <- as.matrix(wide[, -"HugoSymbol"])
rownames(assay) <- wide$HugoSymbol

out <- list(
  assay = assay,
  gene_meta = gene_meta
)

message("Saving preprocessed mutation object to ", OUTPUT$preprocessedMutation)
saveRDS(out, file = OUTPUT$preprocessedMutation)
