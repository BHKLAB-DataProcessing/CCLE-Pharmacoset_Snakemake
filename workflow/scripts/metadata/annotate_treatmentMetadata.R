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
        file.path("resources", paste0(snakemake@rule, ".RData"))
    )
}
load("resources/annotate_treatmentMetadata.RData")

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
# message("\n\nAnnotating with external annotations")
# annotations <- c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code')

# lapply(seq_along(annotations), function(i) {
#     message(paste0("Annotating with ", annotations[i]))
#     treatment_annotations[
#         !is.na(pubchem.CID), 
#         paste0("pubchem.", annotations[i]) := AnnotationGx::annotatePubchemCompound(pubchem.CID, heading = annotations[i])
#         ]
# })

# gsub every column name to remove the " " 
# colnames(treatment_annotations) <- gsub(" ", "_", colnames(treatment_annotations))

# treatmentMetadata <- merge(
#     treatmentMetadata, 
#     treatment_annotations, 
#     by.x = "CCLE.treatmentid", 
#     by.y = "pubchem.name", all.x = TRUE)


unichem_sources <- AnnotationGx::getUnichemSources(T)
data.table::setkey(unichem_sources, Name)

sources_of_interest <- c("chembl", "drugbank", "chebi", "phamgkb", "lincs", "clinicaltrials", "nih_ncc", "fdasrs", "pharmgkb", "rxnorm")

sourceID <- unichem_sources[Name == "pubchem", SourceID]

message("\n\nAnnotating with unichem...")
annotations <- lapply(treatment_annotations$pubchem.CID, function(x){
  tryCatch({
    result <- AnnotationGx::queryUnichemCompound(type = "sourceID", compound = x, sourceID = sourceID)

    subset <- result$External_Mappings[Name %in% sources_of_interest, .(compoundID, Name)]
    # make Name the column names and the values the compoundID 
    subset$cid <- x
    dcast(subset, cid ~ Name, value.var = "compoundID", fun.aggregate = list)
  }, error = function(e) NULL)
  } 
  ) |> data.table::rbindlist(fill = T)
show(annotations)



unichem_mappings <- copy(annotations)
# for each column, if its a list then make it a string with a comma separator
for(col in names(unichem_mappings)){
  if(is.list(unichem_mappings[[col]])){
    unichem_mappings[[col]] <- sapply(unichem_mappings[[col]], function(x) paste(x, collapse = ","))
  }
}
# Rename columns like drugbank to unichem.DrugBank etc, using the unichem_sources "NameLabel" column 
names(unichem_mappings) <- paste("unichem", unichem_sources[names(unichem_mappings), gsub(" ", "_", NameLabel)], sep = ".")

all_annotated_treatmentMetadata <- merge(treatment_annotations, unichem_mappings, by.x = "pubchem.CID", by.y = "unichem.NA", all.x = T)

####################

annotated_treatmentMetadata <- copy(all_annotated_treatmentMetadata)
message("\n\nAnnotating with ChEMBL using Unichem-obtained ChEMBL IDs")
chembl_mechanisms_dt <- annotated_treatmentMetadata[, AnnotationGx::getChemblMechanism(unichem.ChEMBL)]

chembl_cols_of_interest <- c(
        "molecule_chembl_id",  "parent_molecule_chembl_id", "target_chembl_id", "record_id", 
        "mechanism_of_action", "mechanism_comment", "action_type"
    )

annotated_treatmentMetadata <- merge(
    annotated_treatmentMetadata, 
    chembl_mechanisms_dt[, ..chembl_cols_of_interest], 
    by.x = "unichem.ChEMBL",
    by.y = "molecule_chembl_id", 
    all.x = TRUE
    )

data.table::setnames(
    annotated_treatmentMetadata, 
    chembl_cols_of_interest, 
    paste0("chembl.", chembl_cols_of_interest), 
    skip_absent = TRUE)

annotated_treatmentMetadata <- annotated_treatmentMetadata[!duplicated(pubchem.CID),]


final_treatmentMetadata <- merge(
    treatmentMetadata, 
    annotated_treatmentMetadata, 
    by.x = "CCLE.treatmentid", 
    by.y = "pubchem.name", 
    all.x = TRUE

)

# 6.0 Write out cleaned data
# --------------------------
message("\n\nWriting out cleaned data to ", OUTPUT$treatmentMetadata)

data.table::fwrite(
    final_treatmentMetadata, 
    file = OUTPUT$treatmentMetadata,
    quote = FALSE,
    sep = "\t"
)
