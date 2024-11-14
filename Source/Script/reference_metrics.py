from pathlib import Path
import os
import pandas as pd
import sys

def get(list_of_fasta_path:Path,output_dir:str) -> None :

    # Handling path:

    output_path:Path = Path(output_dir)

    reference_info:dict = {}

    # Loop through references FASTA files.

    for fasta in list_of_fasta_path :

        reference_name:str = fasta.stem
        headers:list = []
        sequences:str =""

        with open(fasta,"r") as infile :

            fasta_content:list[str] = infile.read().splitlines()

        for line in fasta_content :

            if line.startswith(">") : 
                headers.append(line)
            
            else :
                sequences+= line
        
        chrom_number:int = len(headers)
        genome_size:int =  len(sequences)

        reference_info[reference_name] = [chrom_number,genome_size]

    # Generating pd.Dataframe:

    reference_dataframe:pd.DataFrame = pd.DataFrame.from_dict(reference_info, orient="index", columns=["Chrom_number", "Genome_size_(bp)"])
    reference_dataframe.reset_index(inplace=True)  # Reset index to turn row names into a column.
    reference_dataframe.rename(columns={"index": "Genome"}, inplace=True)  # Rename the new column to 'Genome'.
    reference_dataframe["Genome_size_(Mb)"] = round(reference_dataframe["Genome_size_(bp)"] / 1000000 , 2) 

    # Converting pd.Dataframe to CSV.

    reference_dataframe.to_csv(output_path, 
                index=False, mode='w', 
                sep=';', encoding='utf-8' )

def main() :
    
    try :
        fasta_dir, metric_file = sys.argv[1], sys.argv[2]

        get(fasta_dir,metric_file)

    except IndexError as err :

        print(err)

if __name__ == "__main__" :

    main()