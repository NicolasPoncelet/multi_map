from Script import get_metrics  
import os

# Defining input directories.
FASTQ_DIR = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data/Fastq"
REFERENCE_DIR = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data/Reference"
ANALYSIS_DIR = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data"

# Get from input directories samples and reference names.

fastq_files = [fastq for fastq in os.listdir(FASTQ_DIR) if fastq.endswith(".fastq.gz")]
reference_files = [fasta for fasta in os.listdir(REFERENCE_DIR) if fasta.endswith((".fasta"))]

SAMPLES = list(set([fastq.split("_R")[0] for fastq in fastq_files]))
REFERENCES = [fasta.split(".")[0] for fasta in reference_files]
print(REFERENCES)

rule all:
    input:
        expand("{analysis_dir}/Metrics/{genome}/{sample}.flagstat", 
            analysis_dir=ANALYSIS_DIR, genome=REFERENCES, sample=SAMPLES)

rule bwa_index:
    input:
        ref = f"{REFERENCE_DIR}/{{genome}}.fasta"
    output:
        amb = f'{REFERENCE_DIR}/{{genome}}.fasta.amb',
        ann = f'{REFERENCE_DIR}/{{genome}}.fasta.ann',
        bwt = f'{REFERENCE_DIR}/{{genome}}.fasta.bwt',
        pac = f'{REFERENCE_DIR}/{{genome}}.fasta.pac',
        sa  = f'{REFERENCE_DIR}/{{genome}}.fasta.sa'
    shell:
        """
        bwa index {input.ref}
        """

rule bwa_map:
    input:
        read_1 = f"{FASTQ_DIR}/{{sample}}_R1.fastq.gz",
        read_2 = f"{FASTQ_DIR}/{{sample}}_R2.fastq.gz",
        genome = f"{REFERENCE_DIR}/{{genome}}.fasta",
        index_files = rules.bwa_index.output,  

    output:
        bam = temp(f"{ANALYSIS_DIR}/Mapped_reads/{{genome}}/{{sample}}.bam"),

    shell:
        """
        bwa mem {input.genome} {input.read_1} {input.read_2} | samtools view -Sb - > {output.bam}
        """

rule samtools_sort:
    input:
        bam = rules.bwa_map.output.bam
    output:
        sorted_bam = f"{ANALYSIS_DIR}/Mapped_reads/{{genome}}/{{sample}}.sorted.bam"

    shell:
        """
        samtools sort {input.bam} -o {output.sorted_bam}
        """

rule samtools_index:
    input:
        sorted_bam = rules.samtools_sort.output.sorted_bam
    output:
        index = f"{ANALYSIS_DIR}/Mapped_reads/{{genome}}/{{sample}}.sorted.bam.bai"

    shell:
        """
        samtools index {input.sorted_bam}
        """

rule samtools_flagstat:
    input:
        sorted_bam = rules.samtools_sort.output.sorted_bam,
        indexed_bam = rules.samtools_index.output.index
    output:
        flagstat = f"{ANALYSIS_DIR}/Metrics/{{genome}}/{{sample}}.flagstat"
    
    shell:
        """
        samtools flagstat {input.sorted_bam} > {output.flagstat}
        """