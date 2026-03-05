import os
from sys import argv
from pathlib import Path
from isolate_parser import Isolate
from isolate_parser.tools import helper_functions as hf

ASSEMBLIES_DIR: Path = Path(argv[1])
OUTPUT_DIR: Path = Path(argv[2])
WUSTL_KEY: Path = argv[3]

isolates: list = [
    os.path.basename(SAMPLE) for SAMPLE in ASSEMBLIES_DIR.iterdir() if SAMPLE.is_dir()
]

hf.output_SLURM_config_TSV_checkm(OUTPUT_DIR, ASSEMBLIES_DIR)
hf.create_checkm_script(OUTPUT_DIR, isolates, "fasta", WUSTL_KEY)
