from sys import argv
from pathlib import Path
from isolate_parser import Isolate
from isolate_parser.tools import helper_functions as hf

READS_DIR: Path = Path(argv[1])
OUTPUT_DIR: Path = Path(argv[2])
WUSTL_KEY = argv[3]

sample_ids: list = hf.create_isolate_list(READS_DIR)
isolates: list = []
for sample in sample_ids:
    isolate = Isolate(sample)
    isolates.append(isolate)

for isolate in isolates:
    isolate.search_for_sequencing_reads(READS_DIR)

hf.output_SLURM_config_TSV_fastp(OUTPUT_DIR, isolates)
hf.create_fastp_script(OUTPUT_DIR, isolates, WUSTL_KEY)
