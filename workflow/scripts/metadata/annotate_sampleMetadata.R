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
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)
sampleMetadata <- sampleMetadata[, .(CCLE.depMapID, CCLE.sampleid, CCLE.name)]
depMapIDS <- sampleMetadata[CCLE.depMapID != "", CCLE.depMapID]
options("mc.cores" = THREADS)

(mapped_cells <- AnnotationGx::mapCell2Accession(
    depMapIDS,
    from = "dr",
))


mapped_sampleMetadata <- merge(
    sampleMetadata, 
    mapped_cells[!is.na(accession),], 
    by.x = "CCLE.depMapID", 
    by.y = "query", 
    all.x = TRUE
)

message("Number of failed mappings: ", sum(is.na(mapped_sampleMetadata$accession)))
print(missing <- mapped_sampleMetadata[is.na(accession), ])

message("Trying again using the sample name")
missing[, c("cellLineName", "accession") := {
    parsed_name <- strsplit(CCLE.sampleid, "_") |> 
        purrr::map_chr(1)
    result <- AnnotationGx::mapCell2Accession(parsed_name)
    list(result$cellLineName, result$accession)
}]

annotated_sampleMetadata <- data.table::rbindlist(
    list(mapped_sampleMetadata, missing), fill = TRUE
)[order(CCLE.sampleid)]

data.table::setnames(
    annotated_sampleMetadata,
    c("cellLineName", "accession"),
    c("cellosaurus.cellLineName", "cellosaurus.accession")
    )


cellosaurus_fields <- c("sy", "ca", "sx", "ag", "di", "derived-from-site", "misspelling")

stopifnot(all(cellosaurus_fields %in% AnnotationGx:::.cellosaurus_fields()))

# (map_fields <- AnnotationGx::mapCell2Accession(
#     cellosaurus_accessions[, cellosaurus.accession],
#     from = "ac",
#     to = cellosaurus_fields,
#     BPPARAM = BPPARAM)
# )

# data.table::setnames(map_fields, cellosaurus_fields, paste0("cellosaurus.", cellosaurus_fields), skip_absent = TRUE)

# cellosaurus_annotations <- 
#     merge(
#         map_fields[, query:= NULL],
#         cellosaurus_accessions,
#         by.x = "query:ac",
#         by.y = "cellosaurus.accession",
#     )
# data.table::setnames(cellosaurus_annotations, "query:ac", "cellosaurus.accession")




# data.table::setkeyv(sampleMetadata, "CCLE.depMapID")

# sampleMetadata <- merge(
#     sampleMetadata, 
#     cellosaurus_annotations, 
#     by.x = "CCLE.depMapID", 
#     by.y = "cellosaurus.depMapID", 
#     all.x = TRUE
# )

# Write Output
# ------------
message("Saving sampleMetadata to: ", OUTPUT$sampleMetadata)
dir.create(dirname(OUTPUT$sampleMetadata), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(
    annotated_sampleMetadata, 
    file = OUTPUT$sampleMetadata, 
    quote = TRUE, 
    sep = "\t", 
    na = "NA", 
    col.names = TRUE
)
