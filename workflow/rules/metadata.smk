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
    script:
        scripts / "metadata/preprocessMetadata.R"
