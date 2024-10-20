import os

# Defining input directories.
FASTQ_DIR = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data/Fastq"
REFERENCE_DIR = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data/Reference"
ANALYSIS_DIR = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data"

# Get from input directories samples and reference names.
fastq_files = [fastq for fastq in os.listdir(FASTQ_DIR) if fastq.endswith(".fastq.gz")]
print(fastq_files)
reference_files = [fasta for fasta in os.listdir(REFERENCE_DIR) if fasta.endswith((".fasta"))]
print(reference_files)

SAMPLES = list(set([fastq.split("_R")[0] for fastq in fastq_files]))
REFERENCES = [fasta.split(".")[0] for fasta in reference_files]

rule all:
    input:
        expand("{analysis_dir}/Sorted_bam/{genome}/{sample}.sorted.bam.bai", 
            analysis_dir=ANALYSIS_DIR, genome=REFERENCES, sample=SAMPLES)

rule bwa_index:
    input:
        genome = "{reference_dir}/{genome}.fasta"  # Updated input path to REFERENCE_DIR
    output:
        genome_index = ["{reference_dir}/{genome}.amb", "{reference_dir}/{genome}.ann", 
                        "{reference_dir}/{genome}.bwt", "{reference_dir}/{genome}.pac", 
                        "{reference_dir}/{genome}.sa"]
    params:
        reference_dir = REFERENCE_DIR,
    shell:
        """
        bwa index {input.genome}
        """

rule bwa_map:
    input:
        fastq1 = f"{FASTQ_DIR}/{{sample}}_R1.fastq.gz",
        fastq2 = f"{FASTQ_DIR}/{{sample}}_R2.fastq.gz",
        genome = f"{REFERENCE_DIR}/{{genome}}.fasta",

    output:
        bam = f"{ANALYSIS_DIR}/Mapped_reads/{{genome}}/{{sample}}.bam",

    shell:
        """
        bwa mem {input.genome} {input.fastq1} {input.fastq2} | samtools view -Sb - > {output.bam}
        """

rule samtools_sort:
    input:
        bam = rules.bwa_map.output.bam
    output:
        sorted_bam = f"{ANALYSIS_DIR}/Sorted_bam/{{genome}}/{{sample}}.sorted.bam"

    shell:
        """
        samtools sort -T {input.bam} -O bam -o {output.sorted_bam}
        """

rule samtools_index:
    input:
        sorted_bam = rules.samtools_sort.output.sorted_bam
    output:
        index = f"{ANALYSIS_DIR}/Sorted_bam/{{genome}}/{{sample}}.sorted.bam.bai"

    shell:
        """
        samtools index {input.sorted_bam}
        """
