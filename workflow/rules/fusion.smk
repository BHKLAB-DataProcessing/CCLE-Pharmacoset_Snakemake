from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

fusion_cfg = config["molecularProfiles"]["fusion"]
filtered_name = filename_from_url(fusion_cfg["filtered"]["url"])
unfiltered_name = filename_from_url(fusion_cfg["unfiltered"]["url"])

conda_env = "../envs/r-bioconductor.yaml"


rule download_Fusion:
    output:
        filtered=rawdata / "fusion" / filtered_name,
        unfiltered=rawdata / "fusion" / unfiltered_name,
    params:
        filtered_url=fusion_cfg["filtered"]["url"],
        unfiltered_url=fusion_cfg["unfiltered"]["url"],
    log:
        logs / "fusion" / "download_Fusion.log",
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.filtered}) $(dirname {log})
        url_f="{params.filtered_url}"; out_f="{output.filtered}"
        if [[ "$url_f" =~ ^https?:// ]]; then
            curl -L "$url_f" -o "$out_f"
        else
            cp "$url_f" "$out_f"
        fi
        url_u="{params.unfiltered_url}"; out_u="{output.unfiltered}"
        if [[ "$url_u" =~ ^https?:// ]]; then
            curl -L "$url_u" -o "$out_u"
        else
            cp "$url_u" "$out_u"
        fi
        '''


rule make_Fusion_SE:
    input:
        filtered=rules.download_Fusion.output.filtered,
        unfiltered=rules.download_Fusion.output.unfiltered,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        fusion_se=procdata / "fusion" / "CCLE_Fusion_SE.RDS",
        feature_metadata=procdata / "fusion" / "CCLE_Fusion_feature_metadata.tsv",
        sample_metadata=procdata / "fusion" / "CCLE_Fusion_sample_metadata.tsv",
    params:
        dataset_config=fusion_cfg,
    log:
        logs / "fusion" / "make_Fusion_SE.log",
    conda:
        conda_env,
    threads:
        1,
    script:
        scripts / "fusion" / "make_Fusion_SE.R"
