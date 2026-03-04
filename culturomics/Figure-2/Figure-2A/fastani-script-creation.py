import os
from sys import argv
from pathlib import Path

CONFIG_DIRECTORY: Path = Path(argv[1])
WUSTL_KEY: str = argv[2]

def create_fastani_script (config_directory: Path, wustl_key: str) -> None:
    """Create fastANI script from each directory containing a config.txt"""
    folder_basename: str = os.path.basename(config_directory)
    with open(f"{config_directory}/fastani.sh", "w") as SH:
        SH.write(
            f"#!/bin/bash\n"
            f"#SBATCH --cpus-per-task=8\n"
            f"#SBATCH --mem=32G\n"
            f"#SBATCH --output=SLURM-outputs/fastani.log\n"
            f"#SBATCH --err=SLURM-outputs/fastani.err\n"
            f"#SBATCH --mail-user={wustl_key}@wustl.edu\n"
            f"#SBATCH --mail-type=ALL\n"
            f"\n"
            f"# Perform all against all ANI calculations of genomes\n"
            f'conda run -n fastani@1.34 fastANI -t "${{SLURM_CPUS_PER_TASK}}" --ql config.txt --rl config.txt -o {config_directory}/{folder_basename}.tsv'
        )

for directory in CONFIG_DIRECTORY.iterdir():
    create_fastani_script(directory, WUSTL_KEY)