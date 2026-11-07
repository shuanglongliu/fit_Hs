"""Broken-symmetry ORCA jobs, one directory per 'u'/'d' configuration.

This directory holds orca_hs.inp (high-spin run), orca_bs.inp (broken-symmetry runs), orca.job and
orca.xyz. The broken-symmetry runs start from the orbitals of the high-spin run, so set up and run
the all-'u' configuration first. Then:

    python main.py setup            # create the directories, write orca.inp / orca.job, link orca.xyz and dft.gbw
    python main.py submit
    python main.py check            # SCF convergence
    python main.py backup           # archive unconverged runs as run_N.*, link dft.gbw to their orbitals
    python main.py restart-inputs   # inputs that read the spin state from dft.gbw (no FlipSpin/FinalMs)
    python main.py resubmit         # resubmit the unconverged jobs
    python main.py magmoms [tol]    # check the Loewdin spin populations against expected_moments
    python main.py energies [tol]   # energies.dat (Ha) for `spinh fit isotropic-exchange`
    python main.py help             # all actions
"""

from spinh.dft import CollinearOrcaJobs

jobs = CollinearOrcaJobs(
    configs=["uuuu"],          # then e.g. "uddd", "uddu", "udud", "uduu", "uudd", "uudu", "uuud"
    spins=[0.5, 0.5, 0.5, 0.5],                 # spin-only |Ms| of each site
    effective_spins=[0.5, 0.5, 2.5, 2.5],       # S of each site in energies.dat
    expected_moments=[1.0, 1.0, 1.0, 1.0],      # |Loewdin spin population| of each site
    n_cpus=32,
    mem_per_core=16000,                         # MB
    magmom_tolerance=0.3,
)

if __name__ == "__main__":
    jobs.run()
