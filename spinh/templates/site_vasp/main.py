"""Noncollinear VASP jobs with one moment direction per site.

This directory holds INCAR (with the fields {mstring}, {rwigs}), POSCAR, KPOINTS, vasp.job and
POTCAR. Then:

    python main.py print          # configurations (local frames) and their directories
    python main.py print-global   # the same in the global frame
    python main.py setup          # create the directories and write the inputs
    python main.py submit
    python main.py check          # convergence;  resubmit: resubmit the unconverged jobs
    python main.py restart-test   # restart plan;  restart: do it (jobs.restart_options)
    python main.py energies       # energies.dat (and exchange_energies.yaml);  dipoles: dipoles.dat
    python main.py help           # all actions

Sites are numbered from 0 here (from 1 in the fit inputs). Directions are (theta, phi) in degrees.
"""

from spinh.dft import SiteJobs, recipes

jobs = SiteJobs(
    magnetic_moments=[7, 7, 4, 4, 4, 4, 4, 4, 4, 4],   # size of the MAGMOM vector of each site
    directions_of="magmom",     # the directions below are of the moments ("magmom") or spins ("spin")
)

# Local frames (optional): {site: z direction (theta, phi), or 3x3 matrix whose rows are the axes};
# frame="local" below means directions in these frames, other sites use the global frame
# jobs.set_frames({0: (42.9, 102.5), 1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]]})

# Collinear configurations
jobs.add(recipes.updown([
    "uuuuuuuuuu", "dduuuuuuuu", "uduuuudddd", "uuuudduudd", "uuudduuddu", "uuududdudu",
    "uuudududud", "uddddduuuu", "uuuuuudddd", "uuudududdd", "uuududdudd", "uuuuududud",
]))

# Single-ion anisotropy: every site along the meridians phi = 0, 180, 90, 270 and the equator
# jobs.add_uniform(recipes.great_circles(n_theta=7, n_phi=12), frame="local")

# Exchange matrix of sites 0 and 3: 9 energy pairs, the other sites along -z; `energies` also
# writes them to exchange_energies.yaml for `spinh fit anisotropic-exchange`
# jobs.add_exchange_pair(0, 3, rest=(180, 0), frame="local")

# Explicit configurations, one (theta, phi) per site
# jobs.add([[(90, 0), (180, 0), (180, 0), (90, 0), ...]], frame="local")

# The same direction for one site in every configuration, or the reverse of it
# jobs.set_direction(0, (180, 0), frame="local")
# jobs.flip(2)

# restart: start unfinished jobs from the WAVECAR of the finished job with the nearest directions
jobs.restart_options = dict(max_angle=45, from_neighbor=True, file_name="WAVECAR", copy_it=False)

if __name__ == "__main__":
    jobs.run()
