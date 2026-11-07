"""Broken-symmetry ORCA jobs for the configurations in data.configs.

    python main.py setup            # create the directories, write orca.inp / orca.job, link orca.xyz and dft.gbw
    python main.py submit
    python main.py check            # SCF convergence
    python main.py backup           # archive unconverged runs as run_N.*, link dft.gbw to their orbitals
    python main.py restart-inputs   # inputs that read the spin state from dft.gbw (no FlipSpin/FinalMs)
    python main.py resubmit         # resubmit the unconverged jobs
    python main.py magmoms [tol]    # check the Loewdin spin populations against expected_magmoms
    python main.py energies [tol]   # energies.dat (Ha) for `spinh isotropic-exchange`
"""

import sys

import data
from spinh.dft.orca import BsDftManager

if __name__ == "__main__":
    BsDftManager(data).run(sys.argv[1:], magmom_tolerance=0.3)
