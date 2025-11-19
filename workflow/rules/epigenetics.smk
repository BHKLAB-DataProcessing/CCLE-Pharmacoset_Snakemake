from pathlib import Path

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

methylation_cfg = config["molecularProfiles"]["epigenetics"]["methylation_rrbs"]
rrbs_datasets = sorted(methylation_cfg.keys())
rrbs_input_paths = {
    dataset: rawdata / "epigenetics" / "methylation" / f"{dataset}.tsv"
    for dataset in rrbs_datasets
}

conda_env = "../envs/r-bioconductor.yaml"
dataset_constraint = "|".join(rrbs_datasets)


rule download_methylation_rrbs:
    output:
        data=rawdata / "epigenetics" / "methylation" / "{dataset}.tsv",
    params:
        url=lambda wildcards: methylation_cfg[wildcards.dataset]["url"],
    log:
        logs / "epigenetics" / "download_rrbs_{dataset}.log",
    wildcard_constraints:
        dataset=dataset_constraint,
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.data}) $(dirname {log})
        {{
          echo "[download] Fetching {params.url}"
          curl -L "{params.url}" | gunzip -c > "{output.data}"
        }} > {log} 2>&1
        """


rule make_Methylation_SE:
    input:
        tss_1kb=rrbs_input_paths["tss_1kb"],
        tss_cpg_clusters=rrbs_input_paths["tss_cpg_clusters"],
        cgi_cpg_clusters=rrbs_input_paths["cgi_cpg_clusters"],
        enhancer_cpg_clusters=rrbs_input_paths["enhancer_cpg_clusters"],
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        methylation_se_list=procdata
        / "epigenetics"
        / "CCLE_RRBS_Methylation_SE_list.RDS",
        tss_1kb_matrix=procdata
        / "epigenetics"
        / "CCLE_RRBS_TSS1kb.tsv",
        tss_cpg_clusters_matrix=procdata
        / "epigenetics"
        / "CCLE_RRBS_TSS_CpG_clusters.tsv",
        cgi_cpg_clusters_matrix=procdata
        / "epigenetics"
        / "CCLE_RRBS_CpG_islands.tsv",
        enhancer_cpg_clusters_matrix=procdata
        / "epigenetics"
        / "CCLE_RRBS_enhancer_CpG_clusters.tsv",
    params:
        dataset_metadata=methylation_cfg,
    log:
        logs / "epigenetics" / "make_Methylation_SE.log",
    conda:
        conda_env,
    threads:
        4
    script:
        scripts / "epigenetics" / "make_Methylation_SE.R"
