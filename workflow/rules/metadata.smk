from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

sample_annotation_name = filename_from_url(config["metadata"]["sampleAnnotation"])
treatment_annotation_name = filename_from_url(config["metadata"]["treatmentAnnotation"])
gencode_name = filename_from_url(config["metadata"]["referenceGenome"]["url"])

annotationGx_docker = config["containers"]["annotationGx"]

################################################################################################
# DOWNLOAD RULES
################################################################################################
rule downloadSampleMetadata:
    output:
        sampleMetadata = metadata / sample_annotation_name
    params:
        url=config["metadata"]["sampleAnnotation"]
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.sampleMetadata})
        curl -L "{params.url}" -o "{output.sampleMetadata}"
        """

rule downloadTreatmentMetadata:
    output:
        treatmentAnnotation = metadata / treatment_annotation_name
    params:
        url=config["metadata"]["treatmentAnnotation"]
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.treatmentAnnotation})
        curl -L "{params.url}" -o "{output.treatmentAnnotation}"
        """

rule preprocessMetadata:
    input:
        sampleAnnotation = metadata / sample_annotation_name,
        treatmentAnnotation = metadata / treatment_annotation_name
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
    output:
        CCLE_GENCODE = metadata / "referenceGenome" / gencode_name
    params:
        url=config["metadata"]["referenceGenome"]["url"]
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.CCLE_GENCODE})
        curl -L "{params.url}" -o "{output.CCLE_GENCODE}"
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
