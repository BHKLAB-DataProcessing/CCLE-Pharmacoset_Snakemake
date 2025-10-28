from pathlib import Path
from urllib.parse import unquote, urlparse

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

rppa = config["molecularProfiles"]["rppa"]["harmonized_rppa"]

matrix_filename = rppa.get("filename")
if not matrix_filename:
    matrix_filename = unquote(Path(urlparse(rppa["url"]).path).name)

conda_env = "../envs/r-bioconductor.yaml"


rule download_RPPA:
    output:
        matrix=rawdata / "proteomics" / matrix_filename,
    params:
        url=rppa["url"],
    log:
        logs / "proteomics" / "download_RPPA.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        curl -L "{params.url}" -o {output.matrix} > {log} 2>&1
        """


rule make_RPPA_SE:
    input:
        matrix=rules.download_RPPA.output.matrix,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        rppa_se=procdata / "proteomics" / "CCLE_RPPA_SE.RDS",
        feature_metadata=procdata / "proteomics" / "CCLE_RPPA_feature_metadata.tsv",
        sample_metadata=procdata / "proteomics" / "CCLE_RPPA_sample_metadata.tsv",
    log:
        logs / "proteomics" / "make_RPPA_SE.log"
    conda:
        conda_env
    script:
        scripts / "proteomics" / "make_RPPA_SE.R"
