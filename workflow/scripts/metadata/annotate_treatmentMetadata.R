## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if (exists("snakemake")) {
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  # setup logger if log file is provided
  if (length(snakemake@log) > 0) {
    sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
  }

  save.image(
    file.path("resources", paste0(snakemake@rule, ".RData"))
  )
}
load("resources/annotate_treatmentMetadata.RData")

library(data.table)
message("Starting annotate_treatmentMetadata.R")

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  tolower(value) %in% c("1", "true", "yes", "y")
}

env_number <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || parsed <= 0) {
    return(default)
  }

  parsed
}

run_rscript_with_timeout <- function(script, args, timeout_seconds, label) {
  output_rds <- tempfile(fileext = ".rds")
  stdout_file <- tempfile()
  stderr_file <- tempfile()
  on.exit(unlink(c(output_rds, stdout_file, stderr_file)), add = TRUE)

  status <- withCallingHandlers(
    system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", script, args, output_rds),
      stdout = stdout_file,
      stderr = stderr_file,
      timeout = timeout_seconds
    ),
    warning = function(w) {
      message(conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (!identical(status, 0L)) {
    stdout <- readLines(stdout_file, warn = FALSE)
    stderr <- readLines(stderr_file, warn = FALSE)
    details <- paste(c(stdout, stderr), collapse = " ")
    stop(
      label,
      " failed with status ",
      status,
      if (nzchar(details)) paste0(": ", details) else "",
      call. = FALSE
    )
  }

  readRDS(output_rds)
}

query_unichem_compound <- function(compound, source_id, timeout_seconds) {
  run_rscript_with_timeout(
    script = file.path(
      "workflow",
      "scripts",
      "metadata",
      "query_unichem_compound.R"
    ),
    args = c(compound, source_id),
    timeout_seconds = timeout_seconds,
    label = paste("UniChem lookup for PubChem CID", compound)
  )
}

get_chembl_mechanism <- function(chembl_ids, timeout_seconds) {
  response_dts <- lapply(chembl_ids, function(chembl_id) {
    tryCatch(
      {
        run_rscript_with_timeout(
          script = file.path(
            "workflow",
            "scripts",
            "metadata",
            "query_chembl_mechanism.R"
          ),
          args = chembl_id,
          timeout_seconds = timeout_seconds,
          label = paste("ChEMBL mechanism lookup for", chembl_id)
        )
      },
      error = function(e) {
        message(
          "ChEMBL mechanism lookup failed for ",
          chembl_id,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  })
  response_dts <- response_dts[lengths(response_dts) > 0]

  all_cols <- c(
    "action_type",
    "binding_site_comment",
    "direct_interaction",
    "disease_efficacy",
    "max_phase",
    "mec_id",
    "mechanism_comment",
    "mechanism_of_action",
    "mechanism_refs",
    "molecular_mechanism",
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "record_id",
    "selectivity_comment",
    "site_id",
    "target_chembl_id",
    "variant_sequence"
  )

  if (!length(response_dts)) {
    return(data.table::as.data.table(stats::setNames(
      replicate(length(all_cols), character(), simplify = FALSE),
      all_cols
    )))
  }

  response_dts <- lapply(response_dts, function(x) {
    for (col in setdiff(all_cols, names(x))) {
      x[, (col) := NA]
    }
    x
  })

  data.table::rbindlist(response_dts, fill = TRUE)
}

# 1.0 Read Input Data
# -------------------
message("\n\nINPUT: ", INPUT)
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)
message(
  "treatmentMetadata: ",
  paste0(capture.output(str(treatmentMetadata)), collapse = "\n")
)

# 2.0 Map treatments to PubChem CID
# ---------------------------------
message(
  "\n\nRunning AnnotationGx::mapCompound2CID on ",
  nrow(treatmentMetadata),
  " treatments"
)
(compound_nameToCIDS <- AnnotationGx::mapCompound2CID(
  treatmentMetadata[, CCLE.treatmentid],
  first = TRUE
))

valid_compound_cids <- compound_nameToCIDS[!is.na(cids)]
missing_cid_names <- compound_nameToCIDS[is.na(cids), unique(name)]
if (length(missing_cid_names)) {
  message(
    "\n\nSkipping PubChem property retrieval for ",
    length(missing_cid_names),
    " treatment(s) without a mapped CID: ",
    paste(missing_cid_names, collapse = ", ")
  )
}


# 3.0 Get PubChem Properties
properties <- c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES')
if (nrow(valid_compound_cids) > 0) {
  message(
    "\n\nGetting the following properties from PubChem: ",
    paste(properties, collapse = " "),
    " for ",
    nrow(valid_compound_cids),
    " compounds"
  )
  pubchemProperties <- AnnotationGx::getPubchemCompound(
    ids = unique(valid_compound_cids$cids),
    from = 'cid',
    to = 'property',
    properties = properties
  )
  pubchemProperties <- data.table::as.data.table(pubchemProperties)
} else {
  message("\n\nNo valid CIDs available; skipping PubChem property retrieval.")
  pubchemProperties <- data.table::data.table(
    CID = integer(),
    Title = character(),
    MolecularFormula = character(),
    InChIKey = character(),
    CanonicalSMILES = character()
  )[0]
}

# 4.0 Merge the annotations
message("\n\nMerging the annotations")
(treatment_annotations <- merge(
  compound_nameToCIDS[, cids := as.character(cids)],
  pubchemProperties[, CID := as.character(CID)],
  by.x = "cids",
  by.y = "CID",
  all.x = TRUE
))
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

sources_of_interest <- c(
  "chembl",
  "drugbank",
  "chebi",
  "phamgkb",
  "lincs",
  "clinicaltrials",
  "nih_ncc",
  "fdasrs",
  "pharmgkb",
  "rxnorm"
)

skip_unichem <- env_flag("CCLE_SKIP_UNICHEM")
unichem_timeout_seconds <- env_number("CCLE_UNICHEM_TIMEOUT_SECONDS", 20)
unichem_sources <- NULL
sourceID <- NA_integer_

if (skip_unichem) {
  message("\n\nSkipping UniChem annotation because CCLE_SKIP_UNICHEM is set.")
  annotations_list <- list()
} else {
  unichem_sources <- tryCatch(
    {
      sources <- AnnotationGx::getUnichemSources(T)
      data.table::setkey(sources, Name)
      sources
    },
    error = function(e) {
      message(
        "\n\nUnable to retrieve UniChem sources; continuing without external IDs: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(unichem_sources)) {
    annotations_list <- list()
  } else {
    sourceID <- unichem_sources[Name == "pubchem", SourceID]
    message(
      "\n\nAnnotating with unichem (timeout: ",
      unichem_timeout_seconds,
      " seconds per PubChem CID)..."
    )
    annotations_list <- lapply(
      unique(stats::na.omit(treatment_annotations$pubchem.CID)),
      function(x) {
        tryCatch(
          {
            result <- query_unichem_compound(
              compound = x,
              source_id = sourceID,
              timeout_seconds = unichem_timeout_seconds
            )
            subset <- result$External_Mappings[
              Name %in% sources_of_interest,
              .(compoundID, Name)
            ]
            subset$cid <- x
            data.table::dcast(
              subset,
              cid ~ Name,
              value.var = "compoundID",
              fun.aggregate = list
            )
          },
          error = function(e) {
            message(
              "UniChem lookup failed for PubChem CID ",
              x,
              ": ",
              conditionMessage(e)
            )
            NULL
          }
        )
      }
    )
  }
}

annotations_list <- annotations_list[lengths(annotations_list) > 0]

if (length(annotations_list)) {
  annotations <- data.table::rbindlist(
    annotations_list,
    fill = TRUE,
    use.names = TRUE
  )
  annotations[, cid := as.character(cid)]
  show(annotations)

  unichem_mappings <- copy(annotations)
  for (col in names(unichem_mappings)) {
    if (is.list(unichem_mappings[[col]])) {
      unichem_mappings[[col]] <- sapply(
        unichem_mappings[[col]],
        function(x) paste(x, collapse = ",")
      )
    }
  }
  names(unichem_mappings) <- paste(
    "unichem",
    unichem_sources[names(unichem_mappings), gsub(" ", "_", NameLabel)],
    sep = "."
  )
} else {
  message("\n\nNo Unichem mappings returned; continuing without external IDs.")
  unichem_mappings <- data.table::data.table(`unichem.NA` = character())
}

all_annotated_treatmentMetadata <- merge(
  treatment_annotations,
  unichem_mappings,
  by.x = "pubchem.CID",
  by.y = "unichem.NA",
  all.x = T
)

####################

annotated_treatmentMetadata <- copy(all_annotated_treatmentMetadata)
message("\n\nAnnotating with ChEMBL using Unichem-obtained ChEMBL IDs")
chembl_ids <- unique(stats::na.omit(annotated_treatmentMetadata$unichem.ChEMBL))
skip_chembl <- env_flag("CCLE_SKIP_CHEMBL")
chembl_timeout_seconds <- env_number("CCLE_CHEMBL_TIMEOUT_SECONDS", 20)
if (skip_chembl) {
  message(
    "Skipping ChEMBL mechanism annotation because CCLE_SKIP_CHEMBL is set."
  )
} else if (length(chembl_ids) > 0) {
  message(
    "Running ChEMBL mechanism lookups with a ",
    chembl_timeout_seconds,
    " second per-request timeout."
  )
  chembl_mechanisms_dt <- get_chembl_mechanism(
    chembl_ids,
    timeout_seconds = chembl_timeout_seconds
  )

  chembl_cols_of_interest <- c(
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "target_chembl_id",
    "record_id",
    "mechanism_of_action",
    "mechanism_comment",
    "action_type"
  )

  if (nrow(chembl_mechanisms_dt) > 0) {
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
      skip_absent = TRUE
    )
  } else {
    message("No ChEMBL mechanism annotations returned.")
  }
} else {
  message("No ChEMBL IDs available; skipping mechanism annotation.")
}

annotated_treatmentMetadata <- annotated_treatmentMetadata[
  !duplicated(pubchem.CID),
]


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
