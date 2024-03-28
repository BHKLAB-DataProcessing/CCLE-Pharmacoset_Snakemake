## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

    save.image(
        file.path("resources/", paste0(snakemake@rule, ".RData"))
    )
}

load("resources/annotate_sampleMetadata.RData")
library(data.table)
# 1.0 Read Input Data
# -------------------

message("\n\nINPUT: ", INPUT)
raw <- data.table::fread(INPUT$sampleMetadata, header = TRUE, sep = "\t", na.strings = NA_character_)
message("Number of samples: ", nrow(raw))

sampleMetadata <- raw

# Set this for mapping 
options("mc.cores" = THREADS)

options("log_level" = "INFO")   # AnnotationGx logging level

# 2.0 Annotate Sample Metadata
# ----------------------------

# Using the depMapID over the cell line name since this ensures the exact match
# that depmap would use
message("Mapping cell line names to Cellosaurus accessions")
options("log_level" = "INFO")
(mapped_cells <- AnnotationGx::mapCell2Accession(
    sampleMetadata[CCLE.depMapID != "", CCLE.depMapID],
    from = "dr",
))

mapped_cells <- unique(mapped_cells[!is.na(accession),])
mapped_cells <- merge(sampleMetadata, mapped_cells, by.x = "CCLE.depMapID", by.y = "query", all.x = FALSE)
message("Successful mappings: ", nrow(mapped_cells))

failed <- sampleMetadata[!CCLE.depMapID %in% mapped_cells$CCLE.depMapID, ]

message("Number of failed mappings: ", nrow(failed))
print(missing <- failed)


# Here we parse the sample id on the _ character and use the first part to map to the cellosaurus
# accession number. Also use fuzzy matching to try and get the correct cell line name
message("Trying again using the sample name")
missing[, c("cellLineName", "accession") := {
    parsed_name <- strsplit(CCLE.sampleid, "_") |> 
        purrr::map_chr(1)
    
    # The only manual change is for the "NCIH2330" cell line since it isnt included in misspellings
    parsed_name[parsed_name == "NCIH2330"] <- "H2330" 
    result <- AnnotationGx::mapCell2Accession(parsed_name, fuzzy = T)
    list(result$cellLineName, result$accession)
}]

# Combine the mapped and missing data
annotated_sampleMetadata <- data.table::rbindlist(
    list(
        mapped_cells, 
        missing
    ), fill = TRUE
)[order(CCLE.sampleid)]

# Rename the columns from cellosaurus
data.table::setnames(
    annotated_sampleMetadata,
    c("cellLineName", "accession"),
    c("cellosaurus.cellLineName", "cellosaurus.accession")
)

# Some of the "CCLE.name" values are empty so we fill them in with the first part of the sample id
annotated_sampleMetadata[CCLE.name == "", CCLE.name:= data.table::tstrsplit(CCLE.sampleid, "_", fixed = TRUE)[1]]

annotated_accessions <- AnnotationGx::annotateCellAccession(
    accessions = annotated_sampleMetadata$cellosaurus.accession,
)
message("Number of annotated accessions: ", nrow(annotated_accessions))

message("Number of unique categories: ")
annotated_accessions[, .N, by = "category"]

message("Number of unique sexOfCell: ")
annotated_accessions[, .N, by = "sexOfCell"]



annotated_accessions[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "; "))]
annotated_accessions[, diseases := sapply(diseases, function(x) paste(x, collapse = "; "))]

annotated_accessions[, c("crossReferences", "hierarchy", "comments") := NULL]

names(annotated_accessions) <- paste0("cellosaurus.", names(annotated_accessions))
annotated_accessions <- unique(annotated_accessions)

final_annotated <- merge(
    annotated_sampleMetadata, 
    annotated_accessions,
    by= c("cellosaurus.accession", "cellosaurus.cellLineName"),
    all.x = TRUE
) |> unique()


# Write Output
# ------------
message("Saving sampleMetadata to: ", OUTPUT$sampleMetadata)
dir.create(dirname(OUTPUT$sampleMetadata), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(
    final_annotated, 
    file = OUTPUT$sampleMetadata, 
    quote = TRUE, 
    sep = "\t", 
    na = "NA", 
    col.names = TRUE
)
