from pathlib import Path
import gzip
import os
import pandas as pd
import sys

def count_reads_number(fastq_dir:str,output_dir) -> None :

    # Handling path:
    output_path = Path(output_dir)
    dir_path:Path = Path(fastq_dir)

    isolate_read_numbers:dict = {}
    fastq_list = [fastq for fastq in os.listdir(dir_path) if fastq.endswith(".fastq.gz")]

    # Loop through references FASTQ files.

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

    # Generating pd.Dataframe:

    fastq_dataframe:pd.DataFrame = pd.DataFrame.from_dict(isolate_read_numbers, orient="index", columns=["Read_1", "Read_2"])
    fastq_dataframe.reset_index(inplace=True)  # Reset index to turn row names into a column.
    fastq_dataframe.rename(columns={"index": "Isolat"}, inplace=True)  # Rename the new column to 'Isolat'.
    fastq_dataframe["Total_reads"] = fastq_dataframe["Read_1"] + fastq_dataframe["Read_2"] 

    # Converting pd.Dataframe to CSV.

    fastq_dataframe.to_csv(output_path, 
                index=False, mode='w', 
                sep=';', encoding='utf-8' )

def main() :
    
    try :
        fastq_dir, metric_file = sys.argv[1], sys.argv[2]

        count_reads_number(fastq_dir,metric_file)

    except IndexError as err :

        print(err)

if __name__ == "__main__" :

    main()