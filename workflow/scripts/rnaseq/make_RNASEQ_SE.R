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
load(file.path("resources/make_RNASEQ_SE.RData"))

## ------------------- LOAD INPUT ------------------- ##
# These are the file names that will be loaded

file_names <- c(
    "transcripts_tpm",
    "genes_tpm",
    "genes_counts",
    "genes_rpkm"
)

# Load the data, and assign it to the variable name, and print the structure
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

# This is the general GENCODE annotation
gencode_GR <- rtracklayer::import(INPUT$GENCODE_Annotation)

# This is the GENCODE annotation provided by the CCLE website
# UNUSED
# ccle_GR <- rtracklayer::import(INPUT$CCLE_GENCODE)

########################################
# Create Metadata list to be used 

# General metadata to be used by all the SummarizedExperiments
# When using this metadata, add the FILL_ME_IN values:
    # datatype = "FILL_ME_IN",
    # filename = "FILL_ME_IN",
    # numSamples = "FILL_ME_IN",
    # numGenes = "FILL_ME_IN"
metadata <- list(
    annotation = "rnaseq",
    class = "RangedSummarizedExperiment",
    data_source = list(
        description = "Read counts, TPM & FPKM-values for all sequenced models including cell lines and organoids.",
        tool = "STAR v2.7.9a",
        url = "https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip"
    ),
    date = Sys.Date()
)
########################################
# genes_tpm
message("Creating genes_tpm SummarizedExperiment")
data.table::setkeyv(genes_tpm, "gene_id")
stopifnot(length(setdiff(genes_tpm[, gene_id], gencode_GR$gene_id)) == 0)

# remove the gene_id and transcript_ids columns and convert the rest to a matrix
tpmMatrix <- as.matrix(
    genes_tpm[, !c("gene_id", "transcript_ids"), with=FALSE],
    rownames=genes_tpm[["gene_id"]]
)

# align the rownames of tpmMatrix with the gene_id of gencode_GR
rowRanges <- gencode_GR[gencode_GR$type == "gene", ]
tpmMatrix <- tpmMatrix[rowRanges$gene_id, ]

# create the SummarizedExperiment object
genes_tpm_SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        exprs = tpmMatrix
    ),
    rowRanges = rowRanges,
    colData = data.table::data.table(
        sampleid = colnames(tpmMatrix),
        # make a column called batchid that is full of NAs
        batchid = rep(NA, ncol(tpmMatrix))
    )
)
message("Adding metadata to genes_tpm:")
(genes_tpm_SE@metadata  <- c(
    metadata,
    list(
        datatype = "genes_tpm",
        filename = basename(INPUT$genes_tpm),
        numSamples = ncol(tpmMatrix),
        numGenes = nrow(tpmMatrix)
    )
))

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

message("Adding metadata to transcripts_tpm:")
(transcripts_tpm_SE@metadata <- c(
    metadata,
    list(
        datatype = "transcripts_tpm",
        filename = basename(INPUT$transcripts_tpm),
        numSamples = ncol(transcripts_tpm_SE),
        numGenes = nrow(transcripts_tpm_SE)
    )
))
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

message("Adding metadata to genes_counts:")
(genes_counts_SE@metadata <- c(
    metadata,
    list(
        datatype = "genes_counts",
        filename = basename(INPUT$genes_counts),
        numSamples = ncol(genes_counts_SE),
        numGenes = nrow(genes_counts_SE)
    )
))

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

message("Adding metadata to genes_rpkm:")
(genes_rpkm_SE@metadata <- c(
    metadata,
    list(
        datatype = "genes_rpkm",
        filename = basename(INPUT$genes_rpkm),
        numSamples = ncol(genes_rpkm_SE),
        numGenes = nrow(genes_rpkm_SE)
    )
))

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

# skipping this metadata step for now
# names(all_metadata) <- file_names

# result <- lapply(names(all_metadata), function(x){
#     exclude_ <- c("sessionInfo", "gene_annotation", "annotation", "data_source", "date", "datatype")
#     subset_metadata <- all_metadata[[x]][!names(all_metadata[[x]]) %in% exclude_]
#     return(subset_metadata)
# })
# names(result) <- file_names

# save_metadata <- c(
#     metadata,
#     list(
#         assay_metadata = result
#     )
# )

# message("Saving Metadata to ", OUTPUT$metadata)
# jsonlite::write_json(save_metadata, OUTPUT$metadata, pretty = TRUE)
