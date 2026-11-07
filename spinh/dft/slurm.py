"""Slurm submission."""

import subprocess
import time


def sbatch(path, script):
    """``sbatch script`` in the directory ``path``."""
    print(path)
    subprocess.run(["sbatch", script], cwd=path)
    time.sleep(1)
