from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

# parse config
treatmentResponse = config["treatmentResponse"]
rawdata_name = filename_from_url(treatmentResponse["rawdata"]["url"])
processed_name = filename_from_url(treatmentResponse["processed"]["url"])

conda_env = "../envs/treatmentResponse.yaml"

rule build_treatmentResponseExperiment:
    input:
        rawdata = rawdata / "treatmentResponse" / rawdata_name,
        processed = rawdata / "treatmentResponse" / processed_name,
        treatmentMetadata = procdata / metadata / "annotations" / "CCLE_treatmentMetadata_annotated.tsv",
        sampleMetadata = procdata / metadata / "annotations" / "CCLE_sampleMetadata_annotated.tsv",
    output:
        tre =  results / "treatmentResponse" / "treatmentResponseExperiment.rds",
        raw =   procdata / "treatmentResponse" / "treatmentResponse_raw.tsv",
        profiles =   procdata / "treatmentResponse" / "treatmentResponse_profiles.tsv",
    log:
        log = logs / "treatmentResponse" / "build_treatmentResponseExperiment.log",
    threads: 
        30
    resources:
        mem_mb = 96000
    conda:
        conda_env
    script:
        scripts / "treatmentResponse/build_treatmentResponseExperiment.R"

rule download_treatmentResponse:
    output:
        rawdata = rawdata / "treatmentResponse" / rawdata_name,
        processed = rawdata / "treatmentResponse" / processed_name,
    params:
        rawdata_url=treatmentResponse["rawdata"]["url"],
        processed_url=treatmentResponse["processed"]["url"],
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.rawdata})
        curl -L "{params.rawdata_url}" -o "{output.rawdata}";
        curl -L "{params.processed_url}" -o "{output.processed}";
        """
