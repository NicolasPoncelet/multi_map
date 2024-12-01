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
        "report":[f"{path_to_analysis_dir}/Metrics/final.csv"],
        "temporary_report":[f"Ressources/flagstat.csv"],
        "reference_report":[f"Ressources/references.csv"],
        "fastq_report":[f"Ressources/fastq.csv"],
        "html_report":[f'{path_to_analysis_dir}/Metrics/report.html']}
    

    # Dynamically build ouput dict given yaml values.

    if config_yaml["mapping"]["map_to_genome"] :

        all_combinations = list( product(references_name,samples_name))
        print(all_combinations)
        all_flagstat_generated:list[str] = [f'{path_to_analysis_dir}/{genome}/{sample}.flagstat' for genome,sample in all_combinations ]
        print(*all_flagstat_generated, sep ='\n')

        all_outputs["flagstat"] = []



get_all_outputs()