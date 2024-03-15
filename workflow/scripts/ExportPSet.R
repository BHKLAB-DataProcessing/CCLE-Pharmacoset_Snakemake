# Script that creates a metadata file for the Data Nutrition Label (DNL) for the GDSC2 dataset

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
        file.path("rdata_files/ExportPSet.RData")
    )
}
load("rdata_files/ExportPSet.RData")
suppressPackageStartupMessages(library(PharmacoGx, quietly = TRUE))

# load pset
message(paste("Loading: ", INPUT$pset, sep = "\n\t"))
pset <- readRDS(INPUT$pset)

export_dir <- OUTPUT$export_dir 

############################################################################
# PARSING THE SLOTS OF THE PHARMACOSET OBJECT
############################################################################
message("\nParsing the slots of the PharmacoSet object")
message(show(pset))
message("Slots:\n\t", paste(slotNames(pset), collapse = ",\n\t"))
# PSet
#   treatment `data.frame`: Treatment metadata
#   sample `data.frame`: Sample metadata
#   molecularProfiles `MultiAssayExperiment`: Molecular profiles
#   treatmentResponse `TreatmentResponseExperiment`: Treatment response
#   annotation `list`: 
#   curation `list`:

############################################################################
# annotation
############################################################################

annotation <- slot(pset, "annotation")
annotation <- annotation[-which(names(annotation) %in% c("call","sessionInfo"))]
path_ <- file.path(export_dir, "annotation.json")
message(paste("Writing: ", path_, sep = "\n\t"))
dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(annotation, path_, pretty = TRUE)

############################################################################
# curation
############################################################################
curation <- slot(pset, "curation")
lapply(names(curation), function(x){
    df <- curation[[x]]
    path_ <- file.path(export_dir, "curation", paste0(x, ".tsv"))
    message(paste("Writing: ", path_, sep = "\n\t"))
    dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
    write.table(
        df,
        file = path_,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE 
    )
    return(TRUE)
}) |> unlist() |> all() |> stopifnot()


############################################################################
# treatmentResponse
############################################################################
treatment <- slot(pset, "treatment")
path_ <- file.path(export_dir, "treatment.tsv")
dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
message(paste("Writing: ", path_, sep = "\n\t"))
write.table(
    treatment,
    file = path_,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)

############################################################################
# sample
############################################################################
sample <- slot(pset, "sample")
path_ <- file.path(export_dir, "sample.tsv")
message(paste("Writing: ", path_, sep = "\n\t"))
dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
write.table(
    sample,
    file = path_,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)

############################################################################
# molecularProfiles
############################################################################
molecularProfiles <- slot(pset, "molecularProfiles")
molecularProfiles <- slot(molecularProfiles, "ExperimentList")
molecularProfiles_metadata <- lapply(names(molecularProfiles), function(x){
    se <- molecularProfiles[[x]]
    metadata_ <- se@metadata
    if(metadata_$annotation == "rnaseq"){
        metadata_ <- metadata_[-which(names(metadata_) == "sessionInfo")]
    }
    return(metadata_)
})
names(molecularProfiles_metadata) <- names(molecularProfiles)
path_ <- file.path(export_dir, "molecularProfiles", "metadata.json")
message(paste("Writing: ", path_, sep = "\n\t"))
dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(molecularProfiles_metadata, path_, pretty = TRUE)


assay_list <- BiocParallel::bplapply(names(molecularProfiles), function(x){
    se <- molecularProfiles[[x]]
    subdir <- se@metadata$annotation
    df <- SummarizedExperiment::assay(se, SummarizedExperiment::assayNames(se))
    path_ <- file.path(export_dir, "molecularProfiles", subdir, paste0(x, ".tsv"))
    message(paste("Writing: ", path_, sep = "\n\t"))
    dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
    write.table(
        df[1:10,],
        file = path_,
        sep = "\t",
        quote = FALSE,
        row.names = TRUE
    )
    return(TRUE)},
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS)
) |> unlist() |> all() |> stopifnot()


############################################################################
# treatmentResponse
############################################################################

treatmentResponse <- slot(pset, "treatmentResponse")

treatmentResponse_metadata <- treatmentResponse@metadata[-which(names(treatmentResponse@metadata) == "sessionInfo")]

path_ <- file.path(export_dir, "treatmentResponse", "metadata.json")
message(paste("Writing: ", path_, sep = "\n\t"))
dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(treatmentResponse_metadata, path_, pretty = TRUE)

tr_assays <- CoreGx::assayNames(treatmentResponse)
lapply(tr_assays, function(x){
    df <- CoreGx::assay(treatmentResponse, x)
    path_ <- file.path(export_dir, "treatmentResponse", paste0(x, ".tsv"))
    message(paste("Writing: ", path_, sep = "\n\t"))
    dir.create(dirname(path_), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(df, path_, sep = "\t", quote = FALSE)
    return(TRUE)
}) |> unlist() |> all() |> stopifnot()


message("\n\nDone!")

# Load the ggplot2 library
library(ggplot2)

# Create a simple data frame
data <- data.frame(
  x = 1:10,
  y = 1:10
)

# Create the plot
my_plot <- ggplot(data, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_line()     # Add lines

# Save the plot to a file (change "plot.png" to your desired filename and extension)
ggsave("results/exports/plot.png", plot = my_plot, width = 6, height = 4, dpi = 300)
