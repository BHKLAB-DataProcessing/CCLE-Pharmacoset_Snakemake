from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

rppa = config["molecularProfiles"]["rppa"]["tcpa_rppa500"]


download_name = filename_from_url(rppa["url"])
compression = "none"
archive_member = rppa.get("archive_member")
matrix_filename = rppa.get("filename")

if download_name.endswith(".zip"):
    compression = "zip"
    if not archive_member:
        archive_member = Path(download_name).with_suffix(".tsv").name
    if not matrix_filename:
        matrix_filename = archive_member
else:
    if not matrix_filename:
        matrix_filename = download_name
    if not archive_member:
        archive_member = matrix_filename

if not matrix_filename:
    raise ValueError("Unable to determine RPPA matrix filename.")

conda_env = "../envs/r-bioconductor.yaml"


rule download_RPPA:
    output:
        matrix=rawdata / "proteomics" / matrix_filename,
    params:
        url=rppa["url"],
        compression=compression,
        archive_member=archive_member,
        download_name=download_name,
    log:
        logs / "proteomics" / "download_RPPA.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.matrix}) $(dirname {log})
        tmp_file="$(mktemp)"
        {{
          echo "[download] Fetching {params.url}"
          echo "[download] Saving to $(basename "$tmp_file") from {params.url}"
          curl -L "{params.url}" -o "$tmp_file"
          case "{params.compression}" in
            zip)
              echo "[extract] Unpacking {params.archive_member}"
              unzip -p "$tmp_file" "{params.archive_member}" > "{output.matrix}"
              rm -f "$tmp_file"
              ;;
            none|"")
              mv "$tmp_file" "{output.matrix}"
              ;;
            *)
              echo "Unsupported compression type: {params.compression}" >&2
              exit 1
              ;;
          esac
        }} > {log} 2>&1
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
