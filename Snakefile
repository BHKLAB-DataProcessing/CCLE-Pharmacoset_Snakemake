configfile: "workflow/config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scriptDir = Path("workflow/scripts")
include: "workflow/rules/metadata.smk"
include: "workflow/rules/downloadGENCODE.smk"
include: "workflow/rules/treatmentResponse.smk"
include: "workflow/rules/rnaseq.smk"

# rules complete:
# make_RNASEQ_SE
# annotate_treatmentMetadata
# build_treatmentResponseExperiment

"""
snakemake \
    --profile workflow/profiles/labserver \
    make_RNASEQ_SE build_treatmentResponseExperiment \
    annotate_treatmentMetadata annotate_sampleMetadata
"""

rule all:
    input:
        results / "CCLE_PSet.RDS",

rule build_MultiAssayExperiment:
    input:
        summarizedExperiment_lists = [
                rules.make_RNASEQ_SE.output.rse_list,
            ],
        metadata_list = [
            rules.make_RNASEQ_SE.output.metadata,
            ],
        sampleMetadata = rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        mae = results / "CCLE_MultiAssayExperiment.RDS"
    script:
        scriptDir / "build_MultiAssayExperiment.R"

rule build_PharmacoSet:
    input:
        multiAssayExperiment = rules.build_MultiAssayExperiment.output.mae,
        treatmentResponseExperiment = rules.build_treatmentResponseExperiment.output.tre,
        treatmentMetadata = rules.annotate_treatmentMetadata.output.treatmentMetadata,
        sampleMetadata = rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        pset = results / "CCLE_PSet.RDS"
    script:
        scriptDir / "build_PharmacoSet.R"


