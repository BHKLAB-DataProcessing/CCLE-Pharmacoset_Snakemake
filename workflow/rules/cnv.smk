from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

cnv = config["molecularProfiles"]["cnv"]
cnv_entry = cnv["copynumber_gene_wgs"]
cnv_gz_name = filename_from_url(cnv_entry["url"])
cnv_path = Path(cnv_gz_name)
cnv_name = (
    cnv_path.with_suffix("").name
    if cnv_path.suffix == ".gz"
    else cnv_path.name
)

conda_env = "../envs/r-bioconductor.yaml"

rule downloadCNV:
    output:
        cnv = rawdata / "cnv" / cnv_name,
    params:
        url=cnv_entry["url"]
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.cnv})
        url="{params.url}"
        if [[ "$url" =~ ^https?:// ]]; then
            if [[ "$url" == *.gz ]]; then
                curl -L "$url" | gunzip -c > "{output.cnv}"
            else
                curl -L "$url" -o "{output.cnv}"
            fi
        else
            if [[ "$url" == *.gz ]]; then
                gunzip -c "$url" > "{output.cnv}"
            else
                cp "$url" "{output.cnv}"
            fi
        fi
        """

rule make_CNV_SE:
    input:
        cnv = rawdata / "cnv" / cnv_name,
        sampleMetadata = procdata / metadata / "preprocessed_sampleMetadata.tsv",
        GENCODE_Annotation =
            gencodeAnnotation(
                dirPath=metadata,
                ref_build=config["metadata"]["referenceGenome"]["build"],
                gencode_ver=config["metadata"]["referenceGenome"]["release"],
                species=config["metadata"]["referenceGenome"]["species"]),
    output:
        CNV_SE = procdata / "cnv" / "CNV_SE.RDS"
    log:
        logs / "cnv" / "make_CNV_SE.log"
    conda:
        conda_env
    threads:
        1
    script:
        scripts / "cnv/make_CNV_SE.R"
