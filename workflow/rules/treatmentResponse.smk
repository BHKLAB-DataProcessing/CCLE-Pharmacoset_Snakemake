from pathlib import Path

configfile: "workflow/config/pipeline.yaml"
conda: "workflow/envs/snakemake.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

# parse config
treatmentResponse = config["treatmentResponse"]


storage HTTP: 
    provider = "http",


rule build_treatmentResponseExperiment:
    input:
        rawdata = rawdata / "treatmentResponse" / "CCLE_NP24.2009_Drug_data_2015.02.24.csv",
        processed = rawdata / "treatmentResponse" / "CCLE_GNF_data_090613.xls",
        treatmentMetadata = procdata / metadata / "preprocessed_treatmentMetadata.tsv",
        sampleMetadata = procdata / metadata / "preprocessed_sampleMetadata.tsv",
    output:
        tre = procdata / "treatmentResponse" / "treatmentResponseExperiment.rds",
        raw = procdata / "treatmentResponse" / "treatmentResponse_raw.tsv",
        profiles = procdata / "treatmentResponse" / "treatmentResponse_profiles.tsv",
    log:
        log = logs / "treatmentResponse" / "build_treatmentResponseExperiment.log",
    threads: 
        30
    resources:
        mem_mb = 96000
    script:
        scripts / "treatmentResponse/build_treatmentResponseExperiment.R"

rule download_treatmentResponse:
    input:
        rawdata = storage.HTTP(treatmentResponse["rawdata"]["url"]),
        processed = storage.HTTP(treatmentResponse["processed"]["url"]),
    output:
        rawdata = rawdata / "treatmentResponse" / "CCLE_NP24.2009_Drug_data_2015.02.24.csv",
        processed = rawdata / "treatmentResponse" / "CCLE_GNF_data_090613.xls",
    shell:
        """
        mv {input.rawdata} {output.rawdata};
        mv {input.processed} {output.processed};
        """