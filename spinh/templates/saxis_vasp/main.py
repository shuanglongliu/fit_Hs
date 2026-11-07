"""Noncollinear VASP jobs over the spin-quantization axis (SAXIS) of the whole system.

This directory holds INCAR (with the field {sstring}), POSCAR, KPOINTS, vasp.job and POTCAR. Run
the job wave_dir first: the others start from its WAVECAR. Then:

    python main.py print          # directions (local frame) and their directories
    python main.py print-global   # the same in the global frame
    python main.py setup          # create the directories and write the inputs
    python main.py submit
    python main.py check          # convergence;  resubmit: resubmit the unconverged jobs
    python main.py restart-test   # jobs without an energy within max_energy of e_ref; restart: resubmit them
    python main.py energies       # energies.dat;  all-energies: for every job directory found
    python main.py help           # all actions
"""

import numpy as np

from spinh.dft import SaxisJobs, recipes

jobs = SaxisJobs(
    wave_dir="0",               # job that supplies WAVECAR ("0" is theta = 0)
    e_ref=-1834.36333292,       # energies further than max_energy (eV) from e_ref count as failed
    max_energy=np.inf,
    de0=1.e-8,                  # converged when the last energy change in `output` is below de0
)

# The local frame: its z direction (theta, phi), e.g. the easy axis found by a first sweep, or a
# 3x3 matrix whose rows are the axes
jobs.set_frame([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

# The meridians phi = 0, 180, 90, 270 and the equator of the local frame
jobs.add_uniform(recipes.great_circles(n_theta=7, n_phi=12), frame="local")

# Other sets, e.g. a fine grid near the local z axis:
# jobs.add_uniform(recipes.meridian(phi=0, theta_min=0, theta_max=6, n=4), frame="local")
# jobs.add_uniform(recipes.parallel(theta=90, phi_min=315, phi_max=345, n=7), frame="local")
# jobs.add_uniform([(0, 0)], frame="global")

if __name__ == "__main__":
    jobs.run()
