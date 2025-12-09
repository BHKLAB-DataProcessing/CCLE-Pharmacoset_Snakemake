from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

cnv_cfg = config["molecularProfiles"]["cnv"]["segments_wgs"]
seg_name = filename_from_url(cnv_cfg["url"])

conda_env = "../envs/r-bioconductor.yaml"


rule download_CNV_segments:
    output:
        seg=rawdata / "cnv" / seg_name,
    params:
        url=cnv_cfg["url"],
    log:
        logs / "cnv" / "download_CNV_segments.log",
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.seg}) $(dirname {log})
        if [[ "{params.url}" =~ ^https?:// ]]; then
          curl -L "{params.url}" -o "{output.seg}"
        else
          cp "{params.url}" "{output.seg}"
        fi
        '''


rule make_CNV_segments_RE:
    input:
        seg=rules.download_CNV_segments.output.seg,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        cnv_segments=procdata / "cnv" / "CNV_segments_RE.RDS",
    log:
        logs / "cnv" / "make_CNV_segments_RE.log",
    conda:
        conda_env,
    threads:
        1,
    script:
        scripts / "cnv" / "make_CNV_segments_RE.R"
