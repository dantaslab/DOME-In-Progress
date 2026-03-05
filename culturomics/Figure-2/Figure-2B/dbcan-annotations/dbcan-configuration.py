from pathlib import Path
import os
from sys import argv
from isolate_parser import Isolate
from isolate_parser.tools import helper_functions as hf

ASSEMBLIES_DIR: Path = Path(argv[1])
OUTPUT_DIR: Path = Path(argv[2])

sample_ids: list = [
    os.path.basename(ASSEMBLY).split(".")[0]
    for ASSEMBLY in ASSEMBLIES_DIR.iterdir()
    if "DOME" in ASSEMBLY.name
]

isolates: list = []
for sample in sample_ids:
    isolate = Isolate(sample)
    isolates.append(isolate)

for isolate in isolates:
    isolate.search_for_genome_assembly(ASSEMBLIES_DIR)

hf.output_SLURM_config_TSV(OUTPUT_DIR, isolates)
