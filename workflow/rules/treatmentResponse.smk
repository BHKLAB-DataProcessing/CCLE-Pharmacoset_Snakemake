from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

# parse config
treatmentResponse = config["treatmentResponse"]

conda_env = "../envs/treatmentResponse.yaml"

rule build_treatmentResponseExperiment:
    input:
        rawdata = rawdata / "treatmentResponse" / "CCLE_NP24.2009_Drug_data_2015.02.24.csv",
        processed = rawdata / "treatmentResponse" / "CCLE_GNF_data_090613.xls",
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
    input:
        rawdata = HTTP.remote(treatmentResponse["rawdata"]["url"]),
        processed = HTTP.remote(treatmentResponse["processed"]["url"]),
    output:
        rawdata = rawdata / "treatmentResponse" / "CCLE_NP24.2009_Drug_data_2015.02.24.csv",
        processed = rawdata / "treatmentResponse" / "CCLE_GNF_data_090613.xls",
    shell:
        """
        mv {input.rawdata} {output.rawdata};
        mv {input.processed} {output.processed};
        """