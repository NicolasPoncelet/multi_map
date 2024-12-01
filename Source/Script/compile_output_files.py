from pathlib import Path
from itertools import product
import yaml

def get_all_outputs() -> dict[str:list[Path]] :

    path_to_yaml:Path = Path('../Config/config.yaml').resolve()

    with open(path_to_yaml, 'r') as infile:
        config_yaml = yaml.load(infile, Loader=yaml.SafeLoader)
    
    # Get inputs path from yaml file to generate output dict for snakefile.

    path_to_fastq:Path = Path(config_yaml["fastq_dir"]).resolve()
    path_to_references:Path = Path(config_yaml["reference_dir"]).resolve()
    path_to_analysis_dir:Path = Path(config_yaml["analysis_dir"]).resolve()

    samples_name:list[str] = [fastq.name.split("_")[0] for fastq in path_to_fastq.glob("*_R1.fastq.gz")]
    references_name:list[str] = [fasta.stem for fasta in path_to_references.glob("*.fasta")]

    # Static output in dict:

    all_outputs:dict[str:[list[Path]]] = {
        "reference_report":[f"Ressources/references.csv"],
        "fastq_report":[f"Ressources/fastq.csv"]    
    }
    # Dynamically build ouput dict given yaml values.

    all_combinations = list( product(references_name,samples_name))

    if config_yaml["mapping"]["map_to_genome"] :

        all_flagstat_generated:list[str] = [f'{path_to_analysis_dir}/Metrics/{genome}/{sample}.flagstat' for genome,sample in all_combinations ]

        all_outputs["flagstat"] = all_flagstat_generated
        all_outputs["html_genomic_report"] = [f'{path_to_analysis_dir}/Metrics/report.html']
        all_outputs["temp_genomic_report"] = [f"Ressources/flagstat.csv"]
        all_outputs["genomic_report"] = [f"{path_to_analysis_dir}/Metrics/final.csv"]

    if not config_yaml["mapping"]["map_to_genome"] :

        all_depth_generated:list[str] = [f'{path_to_analysis_dir}/Metrics/{genome}/{sample}.depth' for genome,sample in all_combinations ]

        all_outputs["depth"] = all_depth_generated
        all_outputs["html_gene_report"] = [f'{path_to_analysis_dir}/Metrics/report.html']
        all_outputs["temp_gene_report"] = [f"Ressources/depth.csv"]
        all_outputs["gene_report"] = [f"{path_to_analysis_dir}/Metrics/final.csv"]


get_all_outputs()