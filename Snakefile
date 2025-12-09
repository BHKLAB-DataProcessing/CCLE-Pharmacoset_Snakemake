from pathlib import Path

configfile: "config/pipeline.yaml"

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
include: "workflow/rules/cnv.smk"
include: "workflow/rules/mutation.smk"
include: "workflow/rules/rppa.smk"
include: "workflow/rules/massspec.smk"
include: "workflow/rules/epigenetics.smk"
include: "workflow/rules/metabolomics.smk"
include: "workflow/rules/mirna.smk"
include: "workflow/rules/exon.smk"
include: "workflow/rules/cnv_segments.smk"

# Until the mtime() issue of gcs storage is resolved for directories, we need to explicitly
# specify all the export results.
# this is tedious since we only know this after we go through the ExportPset Script

exports = [
    "annotation.json",
    "curation/sample.tsv",
    "curation/treatment.tsv",
    "curation/tissue.tsv",
    "treatment.tsv",
    "sample.tsv",
    "molecularProfiles/metadata.json",
    "molecularProfiles/rnaseq/rnaseq.transcripts_tpm.tsv",
    "molecularProfiles/rnaseq/rnaseq.genes_tpm.tsv",
    "molecularProfiles/rnaseq/rnaseq.genes_counts.tsv",
    "molecularProfiles/rnaseq/rnaseq.genes_rpkm.tsv",
    "molecularProfiles/cnv/cnv.genes.tsv",
    "molecularProfiles/mut/mut.genes.tsv",
    "treatmentResponse/metadata.json",
    "treatmentResponse/sensitivity.tsv",
    "treatmentResponse/profiles.tsv"
]
conda_env = "workflow/envs/PharmacoSet.yaml"

rule all:
    input:
        pset=results / "CCLE_PharmacoSet.RDS",
    localrule: True
    # output:
    #     exports = [results / "exports" / export for export in exports],
    # log:
    #     logs / "all.log",
    # conda:
    #     conda_env
    # script:
    #     scriptDir / "ExportPSet.R"


rule build_MultiAssayExperiment:
    input:
        summarizedExperiment_lists=[
            rules.make_RNASEQ_SE.output.rse_list,
            rules.make_CNV_SE.output.CNV_SE,
            rules.make_Mutation_SE.output.processedMutationSE,
            rules.make_RPPA_SE.output.rppa_se,
            rules.make_MassSpec_SE.output.massspec_se,
            rules.make_Methylation_SE.output.methylation_se_list,
            rules.make_Metabolomics_SE.output.metabolomics_se,
            rules.make_miRNA_SE.output.mirna_gct_se,
            rules.make_miRNA_SE.output.mirna_mimat_se,
            rules.make_ExonUsage_SE.output.exon_ratio_se,
            rules.make_ExonUsage_SE.output.exon_denom_se,
            rules.make_CNV_segments_RE.output.cnv_segments,
        ],
        sampleMetadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        mae=results / "CCLE_MultiAssayExperiment.RDS",
    log:
        logs / "build_MultiAssayExperiment.log",
    conda:
        conda_env
    script:
        "workflow/scripts/build_MultiAssayExperiment.R"

rule build_PharmacoSet:
    input:
        multiAssayExperiment=rules.build_MultiAssayExperiment.output.mae,
        treatmentResponseExperiment=rules.build_treatmentResponseExperiment.output.tre,
        treatmentMetadata=rules.annotate_treatmentMetadata.output.treatmentMetadata,
        sampleMetadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        pset=results / "CCLE_PharmacoSet.RDS",
    log:
        logs / "build_PharmacoSet.log",
    conda:
        conda_env
    script:
        scriptDir / "build_PharmacoSet.R"


# # Temporary, separate rule for building the report
# rule build_report:
#     input:
#         pset=results / "CCLE_PSet.RDS",
#     output:
#         results / "report.html",
#     script:
#         scriptDir / "report/report.Rmd"
