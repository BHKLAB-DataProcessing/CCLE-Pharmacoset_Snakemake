from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

massspec_cfg = config["molecularProfiles"]["massspec"]["gygi_ccle"]
conda_env = "../envs/r-bioconductor.yaml"


quant_download_name = filename_from_url(massspec_cfg["url"])
sample_info_download_name = filename_from_url(massspec_cfg["sample_info_url"])


rule download_MassSpec:
    output:
        quant=rawdata / "proteomics" / "massspec" / quant_download_name,
        sample_info=rawdata / "proteomics" / "massspec" / sample_info_download_name,
    params:
        quant_url=massspec_cfg["url"],
        sample_info_url=massspec_cfg["sample_info_url"],
    log:
        logs / "proteomics" / "download_MassSpec.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.quant}) $(dirname {log})
        {{
          echo "[download] Fetching {params.quant_url}"
          curl -L "{params.quant_url}" -o "{output.quant}"
          echo "[download] Fetching {params.sample_info_url}"
          curl -L "{params.sample_info_url}" -o "{output.sample_info}"
        }} > {log} 2>&1
        """


rule make_MassSpec_SE:
    input:
        quant=rules.download_MassSpec.output.quant,
        sample_info=rules.download_MassSpec.output.sample_info,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        massspec_se=procdata / "proteomics" / "CCLE_MassSpec_SE.RDS",
        feature_metadata=procdata / "proteomics" / "CCLE_MassSpec_feature_metadata.tsv",
        sample_metadata=procdata / "proteomics" / "CCLE_MassSpec_sample_metadata.tsv",
    params:
        dataset_release=massspec_cfg.get("release"),
        dataset_config=massspec_cfg,
    log:
        logs / "proteomics" / "make_MassSpec_SE.log"
    conda:
        conda_env
    script:
        scripts / "proteomics" / "make_MassSpec_SE.R"
