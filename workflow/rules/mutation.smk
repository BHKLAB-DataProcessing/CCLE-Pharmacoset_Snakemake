from pathlib import Path

configfile: "workflow/config/pipeline.yaml"
rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

mutation = config["molecularProfiles"]["mutation"]

storage HTTP: 
    provider = "http",

conda_env = "../envs/r-bioconductor.yaml"


rule downloadMutation:
    input:
        oncomapAssay = storage.http(mutation["oncomapAssay"]["url"]),
        oncomap = storage.http(mutation["oncomap"]["url"]),
        hybridCapture = storage.http(mutation["hybridCapture"]["url"]),
    output:
        oncomapAssay=rawdata / "mutation" / "CCLE_Oncomap3_Assays_2012-04-09.csv",
        oncomap=rawdata / "mutation" / "CCLE_Oncomap3_2012-04-09.maf",
        hybridCapture=rawdata / "mutation" / "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",
    shell:
        """
        mv {input.oncomapAssay} {output.oncomapAssay} && \
        mv {input.oncomap} {output.oncomap} && \
        mv {input.hybridCapture} {output.hybridCapture}
        """

###############################################################################
# -- Processing Mutation data -- #
###############################################################################

rule preprocessMutation:
    input:
        oncomapAssay=rawdata / "mutation" / "CCLE_Oncomap3_Assays_2012-04-09.csv",
        oncomap=rawdata / "mutation" / "CCLE_Oncomap3_2012-04-09.maf",
        hybridCapture=rawdata / "mutation" / "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",
    output:
        preprocessedMutation = procdata / "mutation" / "preprocessedMutation.RDS"
    log:
        logs / "mutation" / "preprocessMutation.log"
    conda:
        conda_env
    threads:
        1
    script:
        scripts / "mutation" / "preprocessMutation.R"


rule make_Mutation_SE:
    input:
        preprocessedMutation = procdata / "mutation" / "preprocessedMutation.RDS",
    output:
        processedMutationSE = procdata / "mutation" / "Mutation_SE.RDS"
    conda:
        conda_env
    log:
        logs / "mutation" / "make_Mutation_SE.log"
    threads:
        1 
    script:
        scripts / "mutation" / "make_Mutation_SE.R"