import os
from sys import argv
from pathlib import Path
from isolate_parser import Isolate
from isolate_parser.tools import helper_functions as hf
import pandas as pd

CSV_DIRECTORY: Path = Path(argv[1])
ASSEMBLIES_DIR: Path = Path(argv[2])
OUTPUT_DIR: Path = Path(argv[3])

def taxonomic_CSV_parser(CSV: Path) -> list:
    # List to "hold" Isolate objects within
    isolates: list = []
    # Read in CSV and create Isolate objects
    csv = pd.read_csv(Path(CSV), sep=",").dropna()
    dome_data_dict: dict = dict(zip(csv['user_genome'], csv['classification']))
    # Take zipped dataframe columns and extract information from dict
    for sample_name, taxonomic_classification in dome_data_dict.items():
        dome_isolate: Isolate = Isolate(sample_name)
        dome_isolate.append_isolate_metadata({"classification": taxonomic_classification})
        dome_isolate.search_for_genome_assembly(ASSEMBLIES_DIR)
        isolates.append(dome_isolate)
    return isolates
    
def create_config_file(isolates: list, output_directory: Path, filename: str): 
   output_directory.mkdir(parents=True, exist_ok=True)
   with open(f"{output_directory}/{filename}.txt", "w") as TXT:
        for isolate in isolates:
           TXT.write(
              f"{isolate.genome_assembly.filepath}\n"
            )

for CSV in CSV_DIRECTORY.iterdir():
    isolates: list = taxonomic_CSV_parser(CSV)
    print(isolates[0].metadata['classification'])
    directory_name: list = os.path.basename(isolates[0].metadata["classification"]).split(" ")
    create_config_file(isolates, Path(f'{OUTPUT_DIR}/{directory_name[0]}_{directory_name[1]}'), "config")