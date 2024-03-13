
import os
from pathlib import Path

storage HTTP:
    provider = "http", 

# HOW TO USE
# in a rule, set 
# input:
#   gencodeAnnotation = gencodeAnnotation(dirPath, ref_build, gencode_ver, species="human")
# and this code will handle the rest of the pathing

# IMPORTANT: Gencode v22 and below do not have GRCh37 mappings, any downloaded files will be whatever 
# is available, regardless of how the rule is written

def gencodeAnnotation(dirPath, ref_build, gencode_ver, species="human"):
    return Path(dirPath) / species / f"{ref_build}_v{gencode_ver}" / "annotation.gtf"

def gencodeGenome(dirPath, ref_build, gencode_ver, species="human"):
    return Path(f"{dirPath}/{species}/{ref_build}_v{gencode_ver}/genome.fa")

def gencodeTranscriptome(dirPath, ref_build, gencode_ver, species="human"):
    path = Path(os.path.join(dirPath, species, f'{ref_build}_v{gencode_ver}', "transcriptome.fa"))
    return path


rule downloadAllGencode:
    input:
        expand(
            "{dirPath}/{species}/{ref_build}_v{gencode_ver}/{TYPE}",
            TYPE = ["annotation.gtf", "genome.fa", "transcriptome.fa"],
            dirPath = "references",
            species = "human",
            ref_build = ["GRCh37", "GRCh38"],
            gencode_ver = ["23", "33", "37", "40", "44"],
        ),


######## GENCODE ########
def get_gencode_annotation(ref_build, gencode_release):
    if int(gencode_release) > 22:
        if ref_build == "GRCh37":
            ftp = "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/GRCh37_mapping/gencode.v{gencode_release}lift37.annotation.gtf.gz".format(
                gencode_release = gencode_release.strip()
            )
        else:
            ftp = "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/gencode.v{gencode_release}.annotation.gtf.gz".format(
                gencode_release = gencode_release.strip()
            )
    else:
        # if ref_build == "GRCh37":
            # print("GENCODE release 22 and below do not have GRCh37 mappings")
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/gencode.v{gencode_release}.annotation.gtf.gz"
    return "http://" + ftp


rule getGENCODEannotation:
    output:
        gencode_annotation_file="{refDir}/{species}/{ref_build}_v{gencode_release}/annotation.gtf",
    threads:
        1
    params:
       url = lambda wc: get_gencode_annotation(ref_build = wc.ref_build, gencode_release = wc.gencode_release)
    log:
        "logs/{refDir}/{species}/{ref_build}_v{gencode_release}-annotation_download.log"
    shell:
        """
        wget -O {output.gencode_annotation_file}.gz {params.url} > {log} 2>&1
        gunzip -f {output.gencode_annotation_file}.gz
        """


def get_gencode_genome(ref_build, gencode_release, species="human"):
    if int(gencode_release) > 22:
        if ref_build == "GRCh37":
            ftp_genome = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
        else:
            ftp_genome = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/GRCh38.primary_assembly.genome.fa.gz"
    else:
        if ref_build == "GRCh37":
            print("GENCODE release 22 and below do not have GRCh37 mappings")
        if int(gencode_release) ==19:
            ftp_genome = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz"
        elif (int(gencode_release) == 20) or (int(gencode_release) == 21):
            ftp_genome = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/GRCh38.genome.fa.gz"
        elif int(gencode_release) == 22:
            ftp_genome = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/GRCh38.primary_assembly.genome.fa.gz"

    return "http://" + ftp_genome


rule getGENCODEgenome:
    output:
        gencode_genome_file="{refDir}/{species}/{ref_build}_v{gencode_release}/genome.fa",
    threads:
        1
    params:
        url = lambda wc: get_gencode_genome(ref_build = wc.ref_build, gencode_release = wc.gencode_release)
    log:
        "logs/{refDir}/{species}/{ref_build}_v{gencode_release}-genome_download.log"
    shell:
        """
        wget -O {output.gencode_genome_file}.gz {params.url} > {log} 2>&1
        gunzip -f {output.gencode_genome_file}.gz
        """

def get_gencode_transcriptome(ref_build, gencode_release):
    if int(gencode_release) > 22:
        if ref_build == "GRCh37":
            ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/GRCh37_mapping/gencode.v{gencode_release}lift37.transcripts.fa.gz"
        else:
            ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/gencode.v{gencode_release}.transcripts.fa.gz"
    else:
        if ref_build == "GRCh37":
            print("GENCODE release 22 and below do not have GRCh37 mappings")
        ftp = f"ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_release}/gencode.v{gencode_release}.pc_transcripts.fa.gz"
    return "http://" + ftp


rule getGENCODEtranscriptome:
    output:
        gencode_genome_file="{refDir}/{species}/{ref_build}_v{gencode_release}/transcriptome.fa",
    threads:
        1
    params:
        url = lambda wc: get_gencode_transcriptome(ref_build = wc.ref_build, gencode_release = wc.gencode_release)
    log:
        "logs/{refDir}/{species}/{ref_build}_v{gencode_release}-transcriptome_download.log"
    shell:
        """
        wget -O {output.gencode_genome_file}.gz {params.url} > {log} 2>&1
        gunzip -f {output.gencode_genome_file}.gz
        """
        
