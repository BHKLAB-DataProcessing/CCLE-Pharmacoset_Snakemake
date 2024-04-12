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
load(file.path("resources/make_Mutation_SE.RData"))

message("Loading mutation data")
input <- readRDS(INPUT$preprocessedMutation)
assay <- input$assay
data <- data.table::as.data.table(input$data)

# names(data)
#  [1] "Hugo_Symbol"                   "Entrez_Gene_Id"               
#  [3] "Center"                        "NCBI_Build"                   
#  [5] "Chromosome"                    "Start_position"               
#  [7] "End_position"                  "Strand"                       
#  [9] "Variant_Classification"        "Variant_Type"                 
# [11] "Reference_Allele"              "Tumor_Seq_Allele1"            
# [13] "Tumor_Seq_Allele2"             "dbSNP_RS"                     
# [15] "dbSNP_Val_Status"              "Tumor_Sample_Barcode"         
# [17] "Matched_Norm_Sample_Barcode"   "Match_Norm_Seq_Allele1"       
# [19] "Match_Norm_Seq_Allele2"        "Tumor_Validation_Allele1"     
# [21] "Tumor_Validation_Allele2"      "Match_Norm_Validation_Allele1"
# [23] "Match_Norm_Validation_Allele2" "Verification_Status"          
# [25] "Validation_Status"             "Mutation_Status"              
# [27] "Sequencing_Phase"              "Genome_Change"                
# [29] "Annotation_Transcript"         "Transcript_Strand"            
# [31] "cDNA_Change"                   "Codon_Change"                 
# [33] "Protein_Change"                "Other_Transcripts"            
# [35] "Refseq_mRNA_Id"                "Refseq_prot_Id"               
# [37] "SwissProt_acc_Id"              "SwissProt_entry_Id"           
# [39] "Description"                   "UniProt_AApos"                
# [41] "UniProt_Region"                "UniProt_Site"       
mut <- data.table::as.data.table(data)
tmp <- mut[,.(
    Hugo_Symbol, Transcript_Strand, 
    SwissProt_entry_Id, SwissProt_acc_Id, UniProt_AApos, UniProt_Region, UniProt_Site, Description)] 
# if duplicated Hugo_Symbol, keep the first one without a NA Transcript_Strand or SwissProt_entry_Id
tmp <- tmp[!duplicated(Hugo_Symbol) | (is.na(Transcript_Strand) & is.na(SwissProt_entry_Id)),]

# create variable that drops all mutation specific columns
rrangeData <- data[
    , .(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome, Start_position, End_position, Strand)]
rrangeData <- rrangeData[!duplicated(Hugo_Symbol),]
rrangeData <- merge(rrangeData, tmp, by="Hugo_Symbol", all.x=TRUE)


rowRanges <- GenomicRanges::makeGRangesFromDataFrame(
    df = rrangeData,
    keep.extra.columns=TRUE,
    seqnames.field = "Hugo_Symbol",
    start.field = "Start_position",
    end.field = "End_position",
    strand.field = "Strand",
    na.rm=TRUE
)

message("Creating mutation SummarizedExperiment")
# create a SummarizedExperimentObject
mutation_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        exprs=assay
        ),
    rowRanges = rowRanges,
    colData = data.table::data.table(
        sampleid = colnames(assay),
        # make a column called batchid that is full of NAs
        batchid = rep(NA, ncol(assay))
    )
)

# add metadata
message("Adding metadata to mutation SummarizedExperiment")
(mutation_se@metadata <- list(
    annotation = "mut",
    datatype = "genes",
    class = "RangedSummarizedExperiment",
    filename = "CCLE_Oncomap3_Assays_2012-04-09.csv",
    data_source = snakemake@config$molecularProfiles$mutation$oncomapAssay,
    date = Sys.Date(),
    numSamples = ncol(mutation_se),
    numGenes = nrow(mutation_se)
))


show(mutation_se)

se_list <- list(
    mut.genes = mutation_se
)

message("Saving mutation SE to ", OUTPUT$processedMutationSE)
saveRDS(se_list, file=OUTPUT$processedMutationSE)
