from pathlib import Path
import sys

def get(depth_files:list[Path],csv_report:str, temporary_report:str) -> None :
    """
    Processes depth files to calculate gene coverage metrics and writes results to CSV files.

    Parameters
    ----------
    depth_files : list[Path]
        List of paths to the depth files containing coverage information.
    csv_report : str
        Path to the final CSV report file.
    temporary_report : str
        Path to the temporary CSV report file.

    Notes
    -----
    The function calculates the total size, the size covered (non-zero depth), and the percentage
    coverage for each isolate and gene. The results are written to both the final and temporary
    CSV files.
    """

    # Handling path variable.

    temp_report_path:Path = Path(temporary_report)
    csv_path:Path = Path(csv_report)
    csv_content:str = "\t".join(["Isolat","Gene","Size","Size_covered","Coverage"])

    for depth_report in depth_files :

        # Extract isolate and gene name with Pathlib.

        isolate_name:str = depth_report.stem
        gene_name:str =  depth_report.parent.name

        # Initiating variables.

        current_coverage:int = 0
        total_size:int = 0

        with open (depth_report, "r") as infile :

            content = infile.read().splitlines()
        
        for lines in content :

            _, _, depth = lines.split("\t")

            if depth != '0' :

                current_coverage += 1
            
            total_size += 1
        
        coverage:float = (current_coverage / total_size) * 100 if total_size != 0 else 0

        csv_content += f'\n{isolate_name}\t{gene_name}\t{total_size}\t{current_coverage}\t{coverage:.1f}'
    
    paths = [csv_path, temp_report_path]

    for path in paths:

        path.write_text(csv_content) # method from pathlib module.

def main() :
    """
    Main function to parse command-line arguments and execute the `get` function.

    Command-line Arguments
    -----------------------
    depth_dir : str
        Path to the directory containing depth files.
    metric_file : str
        Path to the final CSV report file.
    temp_file : str
        Path to the temporary CSV report file.

    Raises
    ------
    IndexError
        If the required command-line arguments are not provided.
    """
    
    try :
        depth_dir, metric_file, temp_file = sys.argv[1:4]

        get(depth_dir,metric_file,temp_file)

    except IndexError as err :

        print(err)

if __name__ == "__main__" :

    main()