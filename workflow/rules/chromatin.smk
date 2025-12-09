from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

chrom_cfg = config["molecularProfiles"]["chromatin"]["global_mrm"]
chrom_name = filename_from_url(chrom_cfg["url"])

conda_env = "../envs/r-bioconductor.yaml"


rule download_Chromatin:
    output:
        chromatin = rawdata / "chromatin" / chrom_name,
    params:
        url=chrom_cfg["url"],
    log:
        logs / "chromatin" / "download_Chromatin.log",
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.chromatin}) $(dirname {log})
        url="{params.url}"
        if [[ "$url" =~ ^https?:// ]]; then
            curl -L "$url" -o "{output.chromatin}"
        else
            cp "$url" "{output.chromatin}"
        fi
        '''


rule make_Chromatin_SE:
    input:
        data=rules.download_Chromatin.output.chromatin,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        chromatin_se=procdata / "chromatin" / "CCLE_Chromatin_SE.RDS",
        feature_metadata=procdata / "chromatin" / "CCLE_Chromatin_feature_metadata.tsv",
        sample_metadata=procdata / "chromatin" / "CCLE_Chromatin_sample_metadata.tsv",
    params:
        dataset_config=chrom_cfg,
    log:
        logs / "chromatin" / "make_Chromatin_SE.log",
    conda:
        conda_env,
    threads:
        1,
    script:
        scripts / "chromatin" / "make_Chromatin_SE.R"
