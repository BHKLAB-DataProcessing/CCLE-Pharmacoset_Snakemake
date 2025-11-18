from pathlib import Path
from workflow.utils import filename_from_url

rawdata = Path(config["directories"]["rawdata"])
procdata = Path(config["directories"]["procdata"])
metadata = Path(config["directories"]["metadata"])
results = Path(config["directories"]["results"])
logs = Path(config["directories"]["logs"])
scripts = Path("../scripts")

mutation = config["molecularProfiles"]["mutation"]
somatic_name = filename_from_url(mutation["somatic_mutations"]["url"])


conda_env = "../envs/r-bioconductor.yaml"


rule downloadMutation:
    output:
        somatic=rawdata / "mutation" / somatic_name,
    params:
        somatic_url=mutation["somatic_mutations"]["url"],
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.somatic})
        curl -L "{params.somatic_url}" -o "{output.somatic}"
        """

###############################################################################
# -- Processing Mutation data -- #
###############################################################################

rule preprocessMutation:
    input:
        somatic=rawdata / "mutation" / somatic_name,
        sampleMetadata = procdata / metadata / "preprocessed_sampleMetadata.tsv",
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
        GENCODE_Annotation =
            gencodeAnnotation(
                dirPath=metadata,
                ref_build=config["metadata"]["referenceGenome"]["build"],
                gencode_ver=config["metadata"]["referenceGenome"]["release"],
                species=config["metadata"]["referenceGenome"]["species"]),
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
