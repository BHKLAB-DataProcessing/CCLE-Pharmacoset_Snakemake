from pathlib import Path
configfile: "workflow/config/pipeline.yaml"

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

storage HTTP: 
    provider = "http",

annotationGx_docker = config["containers"]["annotationGx"]

################################################################################################
# DOWNLOAD RULES
################################################################################################
rule downloadSampleMetadata:
    input:
        sampleMetadata = storage.HTTP(config["metadata"]["sampleAnnotation"])
    output:
        sampleMetadata = metadata / "sampleAnnotation.txt"
    shell:
        """
        mv {input.sampleMetadata} {output.sampleMetadata}
        """

rule downloadTreatmentMetadata:
    input:
        treatmentMetadata = storage.HTTP(config["metadata"]["treatmentAnnotation"]),
    output:
        treatmentAnnotation = metadata / "treatmentAnnotation.csv"
    shell:
        """
        mv {input.treatmentMetadata} {output.treatmentAnnotation}
        """

rule preprocessMetadata:
    input:
        sampleAnnotation = metadata / "sampleAnnotation.txt",
        treatmentAnnotation = metadata / "treatmentAnnotation.csv"
    output:
        treatmentMetadata = procdata / metadata / "preprocessed_treatmentMetadata.tsv",
        sampleMetadata = procdata / metadata / "preprocessed_sampleMetadata.tsv",
    threads:
        1
    container: 
        annotationGx_docker
    script:
        scripts / "metadata/preprocessMetadata.R"

rule downloadGenomeFiles:
    input:
        CCLE_GENCODE = 
            storage.HTTP(config["metadata"]["referenceGenome"]["url"]),
    output:
        CCLE_GENCODE = metadata / "referenceGenome" / "gencode.v19.genes.v7_model.patched_contigs.gtf.gz"
    shell:
        """
        mv {input.CCLE_GENCODE} {output.CCLE_GENCODE}
        """

################################################################################################
# TREATMENT METADATA RULES
################################################################################################


rule annotate_treatmentMetadata:
    input:
        treatmentMetadata = procdata / metadata / "preprocessed_treatmentMetadata.tsv",
    output:
        treatmentMetadata = procdata / metadata / "annotations" / "CCLE_treatmentMetadata_annotated.tsv",
    log:
        logs / metadata / "annotate_treatmentMetadata.log"
    container: 
        annotationGx_docker
    threads:
        8
    script:
        scripts / metadata / "annotate_treatmentMetadata.R"


################################################################################################
# SAMPLE METADATA RULES
################################################################################################

rule annotate_sampleMetadata:
    input:
        sampleMetadata = procdata / metadata / "preprocessed_sampleMetadata.tsv",
    output:
        sampleMetadata = procdata / metadata / "annotations" / "CCLE_sampleMetadata_annotated.tsv",
    log:
        logs / metadata / "annotate_sampleMetadata.log"
    container: 
        annotationGx_docker
    threads:
        8
    script:
        scripts / metadata / "annotate_sampleMetadata.R"