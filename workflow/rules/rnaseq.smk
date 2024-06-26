from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

rnaseq = config["molecularProfiles"]["rnaseq"]
conda_env = "../envs/r-bioconductor.yaml"

rule download_RNASEQ:
    input:
        genes_rsem = HTTP.remote(rnaseq["rsem-genes_tpm"]["url"]),
        transcripts_rsem = HTTP.remote(rnaseq["rsem-transcripts_tpm"]["url"]),
        genes_counts = HTTP.remote(rnaseq["genes_counts"]["url"]), 
        genes_rpkm = HTTP.remote(rnaseq["genes_rpkm"]["url"]),
    output:
        genes_tpm=rawdata / "rnaseq" / "CCLE_RNAseq_rsem_genes_tpm_20180929.txt",
        transcripts_tpm=rawdata / "rnaseq" / "CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt",
        genes_counts=rawdata / "rnaseq" / "CCLE_RNAseq_genes_counts_20180929.gct.gz",
        genes_rpkm=rawdata / "rnaseq" / "CCLE_RNAseq_genes_rpkm_20180929.gct.gz",
    log:
        logs / "rnaseq" / "download_RNASEQ.log"
    shell:
        """
        gunzip {input.genes_rsem} -c > {output.genes_tpm} && \
        gunzip {input.transcripts_rsem} -c > {output.transcripts_tpm} && \
        mv {input.genes_rpkm} {output.genes_rpkm} && \
        mv {input.genes_counts} {output.genes_counts} > {log} 2>&1 
        """

rule make_RNASEQ_SE:
    input:
        genes_tpm = rawdata / "rnaseq" / "CCLE_RNAseq_rsem_genes_tpm_20180929.txt",
        transcripts_tpm = rawdata / "rnaseq" / "CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt",
        genes_counts = rawdata / "rnaseq" / "CCLE_RNAseq_genes_counts_20180929.gct.gz",
        genes_rpkm = rawdata / "rnaseq" / "CCLE_RNAseq_genes_rpkm_20180929.gct.gz",
        GENCODE_Annotation = 
            # Function made available from downloadGENCODE.smk
            gencodeAnnotation(
                dirPath=metadata,
                ref_build=config["metadata"]["referenceGenome"]["build"],
                gencode_ver=config["metadata"]["referenceGenome"]["release"],
                species=config["metadata"]["referenceGenome"]["species"]),
        CCLE_GENCODE = metadata / "referenceGenome" / "gencode.v19.genes.v7_model.patched_contigs.gtf.gz"
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