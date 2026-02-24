from sys import argv
from pathlib import Path
from isolate_parser import Isolate
from isolate_parser.tools import helper_functions as hf

TRIMMED_READS_DIR: Path = Path(argv[1])
OUTPUT_DIR: Path = Path(argv[2])
WUSTL_KEY = argv[3]

sample_ids: list = hf.create_isolate_list(TRIMMED_READS_DIR)
isolates: list = []
for sample in sample_ids:
    isolate = Isolate(sample)
    isolates.append(isolate)

for isolate in isolates:
    print(f"Now searching reads for {isolate.sample_name}")
    isolate.search_for_sequencing_reads(TRIMMED_READS_DIR)

hf.output_SLURM_config_TSV_fastp(OUTPUT_DIR, isolates)
print(len(isolates))
hf.create_unicycler_script(OUTPUT_DIR, isolates, WUSTL_KEY)
