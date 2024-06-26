from pathlib import Path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

cnv = config["molecularProfiles"]["cnv"]

conda_env = "../envs/r-bioconductor.yaml"

rule downloadCNV:
    input:
        cnv = HTTP.remote(cnv["copynumber_byGene"]["url"]),
    output:
        cnv = rawdata / "cnv" / "CCLE_copynumber_byGene_2013-12-03.txt",
    shell:
        "gunzip {input.cnv} -c > {output.cnv}"

rule make_CNV_SE:
    input:
        cnv = rawdata / "cnv" / "CCLE_copynumber_byGene_2013-12-03.txt",
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