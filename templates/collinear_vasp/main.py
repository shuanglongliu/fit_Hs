"""Collinear VASP jobs, one directory per 'u'/'d' configuration.

    python main.py setup      # create the directories and write the inputs
    python main.py submit
    python main.py check      # convergence;  resubmit: resubmit the unconverged jobs
    python main.py magmoms    # site moments
    python main.py energies   # energies.dat for `spinh isotropic-exchange`
"""

import sys

import data
from spinh.dft.collinear import CollinearVaspJobs

if __name__ == "__main__":
    # There are 12 exchange couplings. At least 13 magnetic configurations are needed to determine them.
    directories = ["uuuuuuu", "uduuuuu", "uuduuuu", "uuuduuu", "udduuuu", "ududuuu", "uudduuu", "uuuuduu", "uuuuudu", "uuuuuud", "uuuuddu", "uuuudud", "uuuuudd"]
    magmoms = [2, 2, 2, 2, 2, 2, 2]
    effective_spins = [1, 1, 1, 1, 1, 1, 1]

    CollinearVaspJobs(data, directories, magmoms, effective_spins).run(sys.argv[1:])
