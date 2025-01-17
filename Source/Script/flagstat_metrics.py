from pathlib import Path
import sys
from re import findall
import pandas as pd

def get(list_of_flagstat,user_report:str,temporary_report:str) -> None :
    """
    Parse flagstat files and compile the data into a final CSV report.

    Parameters
    ----------
    list_of_flagstat : list of Path
        Paths to individual flagstat report files.
    user_report : str
        Path to the final report CSV file.
    temporary_report : str
        Path to the temporary report CSV file used for the html file.
    """
    
    flagstat_dictionaries:list[dict] = []

    regex_dico:dict = {
            "Total_reads" : r"(\d+ [\+] \d+)( in total \(QC-passed reads \+ QC-failed reads\))",
            "Primary_reads" : r"(\d+ [\+] \d+)( primary[\n])",
            "Secondary_reads" : r"(\d+ [\+] \d+)( secondary[\n])",
            "Supplementary_reads" :r"(\d+ [\+] \d+)( supplementary[\n])",
            "Duplicates_reads" : r"(\d+ [\+] \d+)( duplicates[\n])",
            "Primary_duplicates" : r"(\d+ [\+] \d+)( primary duplicates[\n])",
            "Reads_mapped" :r"(\d+ [\+] \d+)( mapped \(\d+.\d+\% : .+\)[\n])",
            "Reads_primary_mapped": r"(\d+ [\+] \d+)( primary mapped \(\d+.\d+\% : .+\)[\n])",
            "Reads_paired_in_sequencing" :r"(\d+ [\+] \d+)( paired in sequencing[\n])",
            "read1": r"(\d+ [\+] \d+)( read1[\n])",
            "read2" :r"(\d+ [\+] \d+)( read2[\n])",
            "Reads_properly_paired": r"(\d+ [\+] \d+)( properly paired \(\d+.\d+\% : .+\)[\n])",
            "With_itself_and_mate_mapped" :r"(\d+ [\+] \d+)( with itself and mate mapped[\n])",
            "Singletons": r"(\d+ [\+] \d+)( singletons \(\d+.\d+\% : .+\)[\n])",
            "Mate_mapped_to_different_chr" :r"(\d+ [\+] \d+)( with mate mapped to a different chr[\n])"
            }

    path_to_final_report:Path = Path(user_report)
    path_to_temp_report:Path = Path(temporary_report)

    for path_to_flagstat in list_of_flagstat :

        current_dict:dict = {key:0 for key in regex_dico.keys() }

        current_dict["Isolat"] = path_to_flagstat.stem 
        current_dict["Reference"] = path_to_flagstat.parent.name

        with open(path_to_flagstat, 'r') as infile :

            flagstat_content = infile.readlines()

        for keys, patterns in regex_dico.items() :

            for line in flagstat_content :

                if findall(patterns, line) :

                    data_to_catch_plus= findall(patterns, line)[0][0] # Output is a tuple inside a list.
                    data_to_catch, zero = data_to_catch_plus.split("+")

                    if int(zero) != 0:

                        print(f'Warning: {path_to_flagstat.stem} != {zero}')

                    current_dict[keys] = int(data_to_catch)

        flagstat_dictionaries.append(current_dict)
    
    df_flagstat = pd.DataFrame.from_dict(flagstat_dictionaries)
    df_flagstat["Primary_paired_percentage"] = round(df_flagstat['Reads_primary_mapped'] / df_flagstat['Total_reads'] * 100,2)
    df_flagstat["Properly_paired_percentage"] = round(df_flagstat['Reads_properly_paired'] / df_flagstat['Primary_reads'] * 100,2)
    
    df_flagstat.to_csv(path_to_final_report, 
                    index=False, mode='w', 
                    header= True, 
                    sep=';', encoding='utf-8' )
    
    df_flagstat.to_csv(path_to_temp_report, # Generating a second report for HTML report in ressource folder.
                    index=False, mode='w', 
                    header=True, 
                    sep=';', encoding='utf-8' )


if __name__ == "__main__" :
    
    list_of_flagstat, user_report, temporary_report = sys.argv[1], sys.argv[2], sys.argv[3]

    get(list_of_flagstat,user_report,temporary_report)
