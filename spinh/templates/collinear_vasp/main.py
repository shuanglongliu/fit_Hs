"""Collinear VASP jobs, one directory per 'u'/'d' configuration.

This directory holds INCAR (with the fields {magmom}, {nupdown}, {U}, {J}), POSCAR, KPOINTS,
vasp.job and POTCAR. Then:

    python main.py print      # the configurations
    python main.py setup      # create the directories and write the inputs
    python main.py submit
    python main.py check      # convergence;  resubmit: resubmit the unconverged jobs
    python main.py magmoms    # site moments
    python main.py energies   # energies.dat for `spinh fit isotropic-exchange`
    python main.py help       # all actions
"""

from spinh.dft import CollinearVaspJobs

# 12 exchange couplings need at least 13 configurations
jobs = CollinearVaspJobs(
    configs=["uuuuuuu", "uduuuuu", "uuduuuu", "uuuduuu", "udduuuu", "ududuuu", "uudduuu",
             "uuuuduu", "uuuuudu", "uuuuuud", "uuuuddu", "uuuudud", "uuuuudd"],
    magmoms=[2, 2, 2, 2, 2, 2, 2],              # MAGMOM of each site (sign from 'u'/'d')
    effective_spins=[1, 1, 1, 1, 1, 1, 1],      # S of each site in energies.dat
)

if __name__ == "__main__":
    jobs.run()
