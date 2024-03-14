from pathlib import Path

configfile: "workflow/config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

cnv = config["molecularProfiles"]["cnv"]

storage HTTP: 
    provider = "http",

rule downloadCNV:
    input:
        cnv = storage.http(cnv["copynumber_byGene"]["url"]),
    output:
        cnv = rawdata / "cnv" / "CCLE_copynumber_byGene_2013-12-03.txt",
    shell:
        "gunzip {input.cnv} -c > {output.cnv}"


rule make_CNV_SE:
    input:
        cnv= rawdata / "cnv" / "CCLE_copynumber_byGene_2013-12-03.txt",
    output:
        CNV_SE = procdata / "cnv" / "CNV_SE.RDS"
    threads:
        1
    log:
        logs / "cnv" / "make_CNV_SE.log"
    script:
        scripts / "cnv/make_CNV_SE.R"