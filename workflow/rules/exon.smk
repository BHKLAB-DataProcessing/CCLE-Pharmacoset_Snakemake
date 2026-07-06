from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

exon_cfg = config["molecularProfiles"]["exon"]
ratio_name = filename_from_url(exon_cfg["exon_usage_ratio"]["url"])
denom_name = filename_from_url(exon_cfg["exon_usage_denominator"]["url"])

conda_env = "../envs/r-bioconductor.yaml"


rule download_exon_usage:
    output:
        ratio=rawdata / "exon" / ratio_name,
        denom=rawdata / "exon" / denom_name,
    params:
        ratio_url=exon_cfg["exon_usage_ratio"]["url"],
        denom_url=exon_cfg["exon_usage_denominator"]["url"],
    log:
        logs / "exon" / "download_exon_usage.log",
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.ratio}) $(dirname {log})
        {{
          python3 workflow/scripts/download_resource.py --decompress-gzip \
            --source "{params.ratio_url}" \
            --output "{output.ratio}"
          python3 workflow/scripts/download_resource.py --decompress-gzip \
            --source "{params.denom_url}" \
            --output "{output.denom}"
        }} > {log} 2>&1
        '''


rule make_ExonUsage_SE:
    input:
        ratio=rules.download_exon_usage.output.ratio,
        denom=rules.download_exon_usage.output.denom,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        exon_ratio_se=procdata / "exon" / "CCLE_ExonUsageRatio_SE.RDS",
        exon_denom_se=procdata / "exon" / "CCLE_ExonUsageDenom_SE.RDS",
    log:
        logs / "exon" / "make_ExonUsage_SE.log",
    conda:
        conda_env,
    threads:
        1,
    script:
        scripts / "exon" / "make_ExonUsage_SE.R"
