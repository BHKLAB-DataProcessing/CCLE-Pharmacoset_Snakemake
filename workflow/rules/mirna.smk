from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
logs = Path(config["directories"]["logs"])
metadata = Path(config["directories"]["metadata"])
scripts = Path("../scripts")

mirna_cfg = config["molecularProfiles"]["mirna"]
gct_name = filename_from_url(mirna_cfg["mirna_gct"]["url"])
mimat_name = filename_from_url(mirna_cfg["mirna_mimat"]["url"])

conda_env = "../envs/r-bioconductor.yaml"


rule download_miRNA:
    output:
        gct=rawdata / "mirna" / gct_name,
        mimat=rawdata / "mirna" / mimat_name,
    params:
        gct_url=mirna_cfg["mirna_gct"]["url"],
        mimat_url=mirna_cfg["mirna_mimat"]["url"],
    log:
        logs / "mirna" / "download_miRNA.log",
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {output.gct}) $(dirname {log})
        {{
          python3 workflow/scripts/download_resource.py \
            --source "{params.gct_url}" \
            --output "{output.gct}"
          python3 workflow/scripts/download_resource.py \
            --source "{params.mimat_url}" \
            --output "{output.mimat}"
        }} > {log} 2>&1
        '''


rule make_miRNA_SE:
    input:
        gct=rules.download_miRNA.output.gct,
        mimat=rules.download_miRNA.output.mimat,
        sample_metadata=rules.annotate_sampleMetadata.output.sampleMetadata,
    output:
        mirna_gct_se=procdata / "mirna" / "CCLE_miRNA_gct_SE.RDS",
        mirna_mimat_se=procdata / "mirna" / "CCLE_miRNA_mimat_SE.RDS",
    log:
        logs / "mirna" / "make_miRNA_SE.log",
    conda:
        conda_env,
    threads:
        1,
    script:
        scripts / "mirna" / "make_miRNA_SE.R"
