configfile: "workflow/config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])

include: "workflow/rules/metadata.smk"
include: "workflow/rules/downloadGENCODE.smk"
include: "workflow/rules/treatmentResponse.smk"
include: "workflow/rules/rnaseq.smk"

# rules complete:
# make_RNASEQ_SE
# annotate_treatmentMetadata
# build_treatmentResponseExperiment

# snakemake --profile workflow/profiles/labserver  make_RNASEQ_SE build_treatmentResponseExperiment annotate_treatmentMetadata