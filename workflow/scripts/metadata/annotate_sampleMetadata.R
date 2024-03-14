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
        file.path("rdata_files/", paste0(snakemake@rule, ".RData"))
    )
}

load("rdata_files/annotate_sampleMetadata.RData")

# 1.0 Read Input Data
# -------------------

message("\n\nINPUT: ", INPUT)
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)


BPPARAM <- BiocParallel::MulticoreParam(workers = 30, progressbar = TRUE)

(mapped_cells <- AnnotationGx::mapCell2Accession(
    sampleMetadata[!CCLE.depMapID == "", CCLE.depMapID],
    from = "dr",
    BPPARAM= BPPARAM
))

# filter out NA
cellosaurus_accessions <- mapped_cells[!is.na(mapped_cells$ac), .(sampleid = id, accession = ac, depMapID = `query:dr`)]
names(cellosaurus_accessions) <- paste0("cellosaurus.", names(cellosaurus_accessions))

cellosaurus_fields <- c("sy", "ca", "sx", "ag", "di", "derived-from-site", "misspelling")

stopifnot(all(cellosaurus_fields %in% AnnotationGx:::.cellosaurus_fields()))
(map_fields <- AnnotationGx::mapCell2Accession(
    cellosaurus_accessions[, cellosaurus.accession],
    from = "ac",
    to = cellosaurus_fields,
    BPPARAM = BPPARAM)
)

data.table::setnames(map_fields, cellosaurus_fields, paste0("cellosaurus.", cellosaurus_fields), skip_absent = TRUE)

cellosaurus_annotations <- 
    merge(
        map_fields[, query:= NULL],
        cellosaurus_accessions,
        by.x = "query:ac",
        by.y = "cellosaurus.accession",
    )
data.table::setnames(cellosaurus_annotations, "query:ac", "cellosaurus.accession")




data.table::setkeyv(sampleMetadata, "CCLE.depMapID")

sampleMetadata <- merge(
    sampleMetadata, 
    cellosaurus_annotations, 
    by.x = "CCLE.depMapID", 
    by.y = "cellosaurus.depMapID", 
    all.x = TRUE
)

# Write Output
# ------------
message("Saving sampleMetadata to: ", OUTPUT$sampleMetadata)
data.table::fwrite(
    sampleMetadata, 
    file = OUTPUT$sampleMetadata, 
    quote = TRUE, 
    sep = "\t", 
    na = "NA", 
    col.names = TRUE
)
