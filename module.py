from pathlib import Path
import gzip
import os
import pandas as pd

def count_reads_number(fastq_dir:str) -> None :

    # Handling path and file names:

    dir_path:Path = Path(fastq_dir)
    isolate_read_numbers:dict = {}
    fastq_list = [fastq for fastq in os.listdir(dir_path) if fastq.endswith(".fastq.gz")]

    for fastq in fastq_list :

        # Initiating variables:

        read_count:int = 0
        path_to_fastq:Path = dir_path / fastq
        isolate_name:str = path_to_fastq.name.split("_R")[0] # Names are formated like IS001_R1.fastq.gz. Take left part before '_R'.
        pair_read:str = path_to_fastq.name.split("_R")[1][0] # right part after '_R' and first element. either 1 or 2.

        with gzip.open(path_to_fastq,mode="rt") as infile :

            for line in infile :

                if line.startswith("@") :

                    read_count +=1

            if isolate_name not in isolate_read_numbers.keys() :

                isolate_read_numbers[isolate_name] = [None,None]
            
            index_position:int = 0 if pair_read == "1" else 1 # if R1 add in first element of list else second element of list.
            isolate_read_numbers[isolate_name][index_position] = read_count
    
    read_dataframe:pd.DataFrame = pd.DataFrame.from_dict(isolate_read_numbers, orient="index", columns=["Read_1", "Read_2"])
    read_dataframe.reset_index(inplace=True)  # Reset index to turn row names into a column.
    read_dataframe.rename(columns={"index": "Isolat"}, inplace=True)  # Rename the new column to 'Isolat'.
    read_dataframe["Total_reads"] = read_dataframe["Read_1"] + read_dataframe["Read_2"] 

    return read_dataframe

def get_reference_info(fasta_dir:str) -> None :

    # Handling path and file names:

    dir_path:Path = Path(fasta_dir)
    reference_info:dict = {}
    fasta_list = [fastq for fastq in os.listdir(dir_path) if fastq.endswith((".fasta",".fas",".fna",".fa"))]

    for fasta in fasta_list :

        path_to_fasta:Path = dir_path / fasta
        reference_name:str = path_to_fasta.stem
        headers:list = []
        sequences:str =""

        with open(path_to_fasta,"r") as infile :

            fasta_content:list[str] = infile.read().splitlines()

        for line in fasta_content :

            if line.startswith(">") : 
                headers.append(line)
            
            else :
                sequences+= line
        
        chrom_number:int = len(headers)
        genome_size:int =  len(sequences)

        reference_info[reference_name] = [chrom_number,genome_size]

    reference_dataframe:pd.DataFrame = pd.DataFrame.from_dict(reference_info, orient="index", columns=["Chrom_number", "Genome_size_(bp)"])
    reference_dataframe.reset_index(inplace=True)  # Reset index to turn row names into a column.
    reference_dataframe.rename(columns={"index": "Genome"}, inplace=True)  # Rename the new column to 'Genome'.
    reference_dataframe["Genome_size_(Mb)"] = round(reference_dataframe["Genome_size_(bp)"] / 1000000 , 2) 

    return reference_dataframe

def main() :
    
    # isolate_read_numbers:dict = {}
    # dir = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data/Fastq"
    # dir_path = Path(dir)

    # isolate_read_numbers = count_reads_number(dir_path)


    dir = "/home/nponcelet/Documents/03-Script/00_Projet_Perso/02_Bioinfo/33_multi_map/Test_data/Reference"
    reference_data = get_reference_info(dir)

    print(reference_data)


if __name__ == "__main__" :

    main()