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
        file.path("rdata_files", paste0(snakemake@rule, ".RData"))
    )
}
load("rdata_files/annotate_treatmentMetadata.RData")

library(data.table)
message("Starting annotate_treatmentMetadata.R")

# 1.0 Read Input Data
# -------------------
message("\n\nINPUT: ", INPUT)
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)
message("treatmentMetadata: ", paste0(capture.output(str(treatmentMetadata)), collapse = "\n"))

# 2.0 Map treatments to PubChem CID
# ---------------------------------
message("\n\nRunning AnnotationGx::mapCompound2CID on ",nrow(treatmentMetadata), " treatments" )
(compound_nameToCIDS <- AnnotationGx::mapCompound2CID(
    treatmentMetadata[, CCLE.treatmentid],
    first = TRUE
))


# 3.0 Get PubChem Properties
properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES')
message(
    "\n\nGetting the following properties from PubChem: ", 
    paste(properties, collapse= " "), " for ", nrow(compound_nameToCIDS), " compounds")


(pubchemProperties <- 
    compound_nameToCIDS[
        1:nrow(compound_nameToCIDS), 
        AnnotationGx::getPubchemCompound(ids = cids, from = 'cid', to = 'property', properties= properties
    )])

# 4.0 Merge the annotations
message("\n\nMerging the annotations")
(treatment_annotations <- merge(
    compound_nameToCIDS[, cids := as.character(cids)], 
    pubchemProperties[, CID := as.character(CID)], 
    by.x= "cids",  by.y = "CID", all.x = TRUE)
)
data.table::setnames(treatment_annotations, "cids", "CID")
names(treatment_annotations) <- paste0("pubchem.", names(treatment_annotations))

# 5.0 Annotate with external
message("\n\nAnnotating with external annotations")
annotations <- c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code')

lapply(seq_along(annotations), function(i) {
    message(paste0("Annotating with ", annotations[i]))
    treatment_annotations[
        !is.na(pubchem.CID), 
        paste0("pubchem.", annotations[i]) := AnnotationGx::annotatePubchemCompound(pubchem.CID, heading = annotations[i])
        ]
})

# gsub every column name to remove the " " 
colnames(treatment_annotations) <- gsub(" ", "_", colnames(treatment_annotations))

treatmentMetadata <- merge(
    treatmentMetadata, 
    treatment_annotations, 
    by.x = "CCLE.treatmentid", 
    by.y = "pubchem.name", all.x = TRUE)


treatmentMetadata[,pubchem.ChEMBL_ID := data.table::tstrsplit(pubchem.ChEMBL_ID, ";", fixed = TRUE)[1]]

chembl_mechanisms_dt <- treatmentMetadata[, AnnotationGx::getChemblMechanism(pubchem.ChEMBL_ID)]

chembl_cols_of_interest <- c(
        "molecule_chembl_id",  "parent_molecule_chembl_id", "target_chembl_id", "record_id", 
        "mechanism_of_action", "mechanism_comment", "action_type"
    )

treatmentMetadata <- merge(
    treatmentMetadata, 
    chembl_mechanisms_dt[, ..chembl_cols_of_interest], 
    by.x = "pubchem.ChEMBL_ID", 
    by.y = "molecule_chembl_id", all.x = TRUE)

data.table::setnames(
    treatmentMetadata, 
    chembl_cols_of_interest, 
    paste0("chembl.", chembl_cols_of_interest), 
    skip_absent = TRUE)

# 6.0 Write out cleaned data
# --------------------------
message("\n\nWriting out cleaned data to ", OUTPUT$treatmentMetadata)

data.table::fwrite(
    treatmentMetadata, 
    file = OUTPUT$treatmentMetadata,
    quote = FALSE,
    sep = "\t"
)
