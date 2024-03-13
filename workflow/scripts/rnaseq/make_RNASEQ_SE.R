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
load(file.path("rdata_files/make_RNASEQ_SE.RData"))

file_names <- c(
    "transcripts_tpm",
    "genes_tpm",
    "genes_counts",
    "genes_rpkm"
)

for (i in seq_along(file_names)) {
    message(paste0("Loading ", file_names[i]))
    assign(file_names[i], data.table::fread(
        INPUT[[file_names[i]]],
        header = TRUE,
        stringsAsFactors = FALSE,
        sep = "\t"
    ))
    str(get(file_names[i]))
}

message("Loading GENCODE")
gencode_GR <- rtracklayer::import(INPUT$GENCODE_Annotation)

ccle_GR <- rtracklayer::import(INPUT$CCLE_GENCODE)

########################################
# Create Metadata list to be used 
metadata <- list(
    data_source = snakemake@config$molecularProfiles$rnaseq,
    annotation = "rnaseq",
    gene_annotation = lapply(snakemake@config$metadata$referenceGenome, as.character),
    date = Sys.Date()
)


########################################
# genes_tpm
message("Creating genes_tpm SummarizedExperiment")
data.table::setkeyv(genes_tpm, "gene_id")
stopifnot(length(setdiff(genes_tpm[, gene_id], gencode_GR$gene_id)) == 0)
tpmMatrix <- as.matrix(
    genes_tpm[, !c("gene_id", "transcript_ids"), with=FALSE],
    rownames=genes_tpm[["gene_id"]]
)
rowRanges <- gencode_GR[gencode_GR$type == "gene", ]

# align the rownames of tpmMatrix with the gene_id of gencode_GR
tpmMatrix <- tpmMatrix[rowRanges$gene_id, ]

genes_tpm_SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        exprs = tpmMatrix
    ),
    rowRanges = rowRanges,
    colData = data.table::data.table(
        sampleid = colnames(tpmMatrix),
        # make a column called batchid that is full of NAs
        batchid = rep(NA, ncol(tpmMatrix))
    ),
)
genes_tpm_SE@metadata <- c(
    metadata,
    list(
        datatype = "genes_tpm",
        filename = basename(INPUT$genes_tpm),
        samples = ncol(tpmMatrix),
        genes = nrow(tpmMatrix),
        sessionInfo = sessionInfo()
    )
)
show(genes_tpm_SE)
rm(tpmMatrix, rowRanges)

########################################
# transcripts_tpm
message("Creating transcripts_tpm SummarizedExperiment")
data.table::setkeyv(transcripts_tpm, "transcript_id")
stopifnot(length(setdiff(transcripts_tpm[, transcript_id], gencode_GR$transcript_id)) == 0)

tpmMatrix <- as.matrix(
    transcripts_tpm[, !c("transcript_id", "gene_id"), with=FALSE],
    rownames = transcripts_tpm[["transcript_id"]]
)

rowRanges <- gencode_GR[gencode_GR$type == "transcript", ]

# align the rownames of tpmMatrix with the transcript_id of gencode_GR
tpmMatrix <- tpmMatrix[rowRanges$transcript_id, ]

transcripts_tpm_SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        exprs = tpmMatrix
    ),
    rowRanges = rowRanges,
    colData = data.table::data.table(
        sampleid = colnames(tpmMatrix),
        # make a column called batchid that is full of NAs
        batchid = rep(NA, ncol(tpmMatrix))
    ),
)
transcripts_tpm_SE@metadata <- c(
    metadata,
    list(
        datatype = "transcripts_tpm",
        filename = basename(INPUT$transcripts_tpm),
        samples = ncol(tpmMatrix),
        genes = nrow(tpmMatrix),
        sessionInfo = sessionInfo()
    )
)
show(transcripts_tpm_SE)
rm(tpmMatrix, rowRanges)

########################################
# genes_counts
message("Creating genes_counts SummarizedExperiment")
data.table::setkeyv(genes_counts, "Name")       # this table calls the gene_id "Name"
stopifnot(length(setdiff(genes_counts[, Name], gencode_GR$gene_id)) == 0)

countsMatrix <- as.matrix(
    genes_counts[, !c("Name", "Description"), with=FALSE],
    rownames = genes_counts[["Name"]]
)

rowRanges <- gencode_GR[gencode_GR$type == "gene" & gencode_GR$gene_id %in% rownames(countsMatrix), ]

# align the rownames of countsMatrix with the gene_id of gencode_GR
countsMatrix <- countsMatrix[rowRanges$gene_id, ]

genes_counts_SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        exprs = countsMatrix
    ),
    rowRanges = rowRanges,
    colData = data.table::data.table(
        sampleid = colnames(countsMatrix),
        # make a column called batchid that is full of NAs
        batchid = rep(NA, ncol(countsMatrix))
    ),
)

genes_counts_SE@metadata <- c(
    metadata,
    list(
        datatype = "genes_counts",
        filename = basename(INPUT$genes_counts),
        samples = ncol(countsMatrix),
        genes = nrow(countsMatrix),
        sessionInfo = sessionInfo()
    )
)
show(genes_counts_SE)
rm(countsMatrix, rowRanges)

########################################
# genes_rpkm
message("Creating genes_rpkm SummarizedExperiment")
data.table::setkeyv(genes_rpkm, "Name")
stopifnot(length(setdiff(genes_rpkm[, Name], gencode_GR$gene_id)) == 0)

rpkmMatrix <- as.matrix(
    genes_rpkm[, !c("Name", "Description"), with=FALSE],
    rownames = genes_rpkm[["Name"]]
)

rowRanges <- gencode_GR[gencode_GR$type == "gene" & gencode_GR$gene_id %in% rownames(rpkmMatrix), ]

# align the rownames of rpkmMatrix with the gene_id of gencode_GR
rpkmMatrix <- rpkmMatrix[rowRanges$gene_id, ]

genes_rpkm_SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        exprs = rpkmMatrix
    ),
    rowRanges = rowRanges,
    colData = data.table::data.table(
        sampleid = colnames(rpkmMatrix),
        # make a column called batchid that is full of NAs
        batchid = rep(NA, ncol(rpkmMatrix))
    ),
)

genes_rpkm_SE@metadata <- c(
    metadata,
    list(
        datatype = "genes_rpkm",
        filename = basename(INPUT$genes_rpkm),
        samples = ncol(rpkmMatrix),
        genes = nrow(rpkmMatrix),
        sessionInfo = sessionInfo()
    )
)

show(genes_rpkm_SE)
rm(rpkmMatrix, rowRanges)

########################################
# SAVE OUTPUT

rse_list <- list( genes_tpm_SE, transcripts_tpm_SE, genes_counts_SE, genes_rpkm_SE)
names(rse_list) <- paste0("rnaseq.", file_names)


message("Saving RSE List to ", OUTPUT$rse_list)
saveRDS(rse_list, file = OUTPUT$rse_list)

file_names <- c(
    "transcripts_tpm",
    "genes_tpm",
    "genes_counts",
    "genes_rpkm"
)

# Save each assay to output files and collect their metadata
all_metadata <- lapply(file_names, function(x){

    message(paste0("Writing ", x, " to ", OUTPUT[[x]]))
    write.table(
        x = SummarizedExperiment::assay(rse_list[[paste0("rnaseq.", x)]]),
        file = OUTPUT[[x]],
        quote = FALSE,
        sep = "\t",
        row.names = TRUE
    )
    return(rse_list[[paste0("rnaseq.", x)]]@metadata)   
})
names(all_metadata) <- file_names

result <- lapply(names(all_metadata), function(x){
    exclude_ <- c("sessionInfo", "gene_annotation", "annotation", "data_source", "date", "datatype")
    subset_metadata <- all_metadata[[x]][!names(all_metadata[[x]]) %in% exclude_]
    return(subset_metadata)
})
names(result) <- file_names

save_metadata <- c(
    metadata,
    list(
        assay_metadata = result
    )
)

message("Saving Metadata to ", OUTPUT$metadata)
jsonlite::write_json(save_metadata, OUTPUT$metadata, pretty = TRUE)
