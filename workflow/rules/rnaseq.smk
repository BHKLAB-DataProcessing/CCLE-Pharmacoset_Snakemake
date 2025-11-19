from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

rnaseq = config["molecularProfiles"]["rnaseq"]
genes_tpm_gz = filename_from_url(rnaseq["rsem-genes_tpm"]["url"])
transcripts_tpm_gz = filename_from_url(rnaseq["rsem-transcripts_tpm"]["url"])
genes_counts_name = filename_from_url(rnaseq["genes_counts"]["url"])
genes_rpkm_name = filename_from_url(rnaseq["genes_rpkm"]["url"])

genes_tpm_name = Path(genes_tpm_gz).with_suffix("").name
transcripts_tpm_name = Path(transcripts_tpm_gz).with_suffix("").name
gencode_name = filename_from_url(config["metadata"]["referenceGenome"]["url"])

conda_env = "../envs/r-bioconductor.yaml"

rule download_RNASEQ:
    output:
        genes_tpm=rawdata / "rnaseq" / genes_tpm_name,
        transcripts_tpm=rawdata / "rnaseq" / transcripts_tpm_name,
        genes_counts=rawdata / "rnaseq" / genes_counts_name,
        genes_rpkm=rawdata / "rnaseq" / genes_rpkm_name,
    params:
        genes_rsem_url=rnaseq["rsem-genes_tpm"]["url"],
        transcripts_rsem_url=rnaseq["rsem-transcripts_tpm"]["url"],
        genes_counts_url=rnaseq["genes_counts"]["url"],
        genes_rpkm_url=rnaseq["genes_rpkm"]["url"],
    log:
        logs / "rnaseq" / "download_RNASEQ.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.genes_tpm}) $(dirname {log})
        {{
          echo "[download] Fetching genes TPM {params.genes_rsem_url}"
          curl -L "{params.genes_rsem_url}" | gunzip -c > "{output.genes_tpm}"
          echo "[download] Fetching transcripts TPM {params.transcripts_rsem_url}"
          curl -L "{params.transcripts_rsem_url}" | gunzip -c > "{output.transcripts_tpm}"
          echo "[download] Fetching genes RPKM {params.genes_rpkm_url}"
          curl -L "{params.genes_rpkm_url}" -o "{output.genes_rpkm}"
          echo "[download] Fetching genes counts {params.genes_counts_url}"
          curl -L "{params.genes_counts_url}" -o "{output.genes_counts}"
        }} > {log} 2>&1
        """

rule make_RNASEQ_SE:
    input:
        genes_tpm = rawdata / "rnaseq" / genes_tpm_name,
        transcripts_tpm = rawdata / "rnaseq" / transcripts_tpm_name,
        genes_counts = rawdata / "rnaseq" / genes_counts_name,
        genes_rpkm = rawdata / "rnaseq" / genes_rpkm_name,
        GENCODE_Annotation = 
            # Function made available from downloadGENCODE.smk
            gencodeAnnotation(
                dirPath=metadata,
                ref_build=config["metadata"]["referenceGenome"]["build"],
                gencode_ver=config["metadata"]["referenceGenome"]["release"],
                species=config["metadata"]["referenceGenome"]["species"]),
        CCLE_GENCODE = metadata / "referenceGenome" / gencode_name
    output:
        rse_list = procdata / "rnaseq" / "CCLE_RNASEQ_rse_list.RDS",
        genes_tpm = procdata / "rnaseq" / "CCLE_RNAseq_genes_tpm.tsv",
        transcripts_tpm = procdata / "rnaseq" / "CCLE_RNAseq_transcripts_tpm.tsv",
        genes_counts = procdata / "rnaseq" / "CCLE_RNAseq_genes_counts.tsv",
        genes_rpkm = procdata / "rnaseq" / "CCLE_RNAseq_genes_rpkm.tsv"
    log:
        logs / "rnaseq" / "make_RNASEQ_SE.log"
    conda:
        conda_env
    threads:
        3
    script:
        scripts / "rnaseq" / "make_RNASEQ_SE.R"
