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
load("rdata_files/build_MultiAssayExperiment.RData")


# Read in metadata
# ----------------
message(paste("Loading: ", INPUT$sampleMetadata, sep = "\n\t"))
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)

# Read in the summarized experiments
# ----------------------------------
message(paste("Loading: ", INPUT$summarizedExperiment_lists, sep = "\n\t"))
se_list <- lapply(INPUT$summarizedExperiment_lists, function(x){
    message(paste("Loading: ", x, sep = "\n\t"))
    readRDS(x) 
    }
) |> unlist()

se_sampleNames <- lapply(
    se_list, SummarizedExperiment::colnames
    ) |> 
    unlist() |> 
    unique()
    
if(!all(se_sampleNames %in% sampleMetadata$CCLE.sampleid)){
    print("Not all sample names are in the sample metadata")

    message("Removing the following sample names from the summarized experiments:\n\t")
    missing <- setdiff(se_sampleNames, sampleMetadata$CCLE.sampleid)
    message("\t",paste(missing, collapse = "\n\t"))

    se_list <- lapply(se_list, function(x){
        x[,!colnames(x) %in% missing]
    })
}


data.table::setkeyv(sampleMetadata, "CCLE.sampleid")

sampleMetadata <- sampleMetadata[se_sampleNames,] |> unique()

# first_col <- c("CCLE.sampleid", "CCLE.name", "cellosaurus.sampleid")
# # order the sample metadata columns 

# Build MultiAssayExperiment
# --------------------------

# create a data frame for coldata including sampleids and batch ids
colData <- as.data.frame(sampleMetadata, row.names = sampleMetadata$CCLE.sampleid)

colData$sampleid <- colData$CCLE.sampleid
colData$batchid <- 1

message(sprintf("Column data has %d rows and %d columns", nrow(colData), ncol(colData)))
str(colData)


summarizedExperimentLists <- data.table::copy(se_list)
summarizedExperimentLists <- sapply(summarizedExperimentLists, function(x){
    x@colData <- MultiAssayExperiment::DataFrame(
        sampleid = colnames(x),
        batchid = rep(NA, ncol(x)),
        row.names = colnames(x)
    )
    x
})
ExpList <- MultiAssayExperiment::ExperimentList(summarizedExperimentLists)
message(paste("ExperimentList:", capture.output(show(ExpList)), sep = "\n\t"))


# Create a sample map for each experiment in the ExperimentList
sampleMapList <- lapply(summarizedExperimentLists, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})
names(sampleMapList) <- names(ExpList)
message(paste("Sample map list:", capture.output(str(sampleMapList)), sep = "\n\t"))

# Metadata List
# go through each experiment, extract the metadata and add it to a list
metadata_list <- lapply(summarizedExperimentLists, function(se){
    metadata_ <- slot(se, "metadata")
    if(metadata_$annotation == "rnaseq"){
        metadata_ <- metadata_[-which(names(metadata_) == "sessionInfo")]
    }
    metadata_
})



mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList),
    metadata = metadata_list
)


# Write Output
# ------------
message("Saving MultiAssayExperiment to: ", OUTPUT$mae)
saveRDS(mae, file = OUTPUT$mae)
