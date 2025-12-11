from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

metabolomics_cfg = config["molecularProfiles"]["metabolomics"]["ccle_2019"]
metabolomics_filename = filename_from_url(metabolomics_cfg["url"])

conda_env = "../envs/r-bioconductor.yaml"


rule download_Metabolomics:
    output:
        data=rawdata / "metabolomics" / metabolomics_filename,
    params:
        url=metabolomics_cfg["url"],
    log:
        logs / "metabolomics" / "download_Metabolomics.log",
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.data}) $(dirname {log})
        {{
          echo "[download] Fetching {params.url}"
          curl -L "{params.url}" -o "{output.data}"
        }} > {log} 2>&1
        """


rule make_Metabolomics_SE:
    input:
        data=rules.download_Metabolomics.output.data,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        metabolomics_se=procdata / "metabolomics" / "CCLE_Metabolomics_SE.RDS",
        feature_metadata=procdata / "metabolomics" / "CCLE_Metabolomics_feature_metadata.tsv",
        sample_metadata=procdata / "metabolomics" / "CCLE_Metabolomics_sample_metadata.tsv",
    params:
        dataset_config=metabolomics_cfg,
    log:
        logs / "metabolomics" / "make_Metabolomics_SE.log",
    conda:
        conda_env
    threads:
        2
    script:
        scripts / "metabolomics" / "make_Metabolomics_SE.R"
