# DFT jobs

`spinh new <template> <dir>` creates a job directory. It holds the inputs shared by all jobs as
plain files, and `main.py`, which lists the spin configurations and runs an action:
`python main.py <action>`. `python main.py help` lists the actions. The job directories, one per
configuration, are created next to `main.py`.

| Template | Spins | Input files (fields filled per job) | Directory names |
|---|---|---|---|
| `collinear_vasp` | up/down per site (`MAGMOM`, `NUPDOWN`) | `INCAR` (`{magmom}`, `{nupdown}`, `{U}`, `{J}`), `POSCAR`, `KPOINTS`, `vasp.job`, `POTCAR` | `uudu…` |
| `collinear_orca` | broken-symmetry ORCA (`FlipSpin`, `FinalMs`) | `orca_hs.inp`, `orca_bs.inp` (see `spinh.dft.orca`), `orca.job`, `orca.xyz` | `uudu…` |
| `saxis_vasp` | one axis for the whole system (`SAXIS`) | `INCAR` (`{sstring}`), `POSCAR`, `KPOINTS`, `vasp.job`, `POTCAR` | `6476401` |
| `site_vasp` | one moment direction per site (`MAGMOM`, `M_CONSTR`, `RWIGS`) | `INCAR` (`{mstring}`, `{rwigs}`), `POSCAR`, `KPOINTS`, `vasp.job`, `POTCAR` | `4319101_0_0` |

POTCAR is linked into the job directories, the other files are copied (INCAR with its fields
filled in). Commented-out lines of INCAR are kept as they are.

The job sets are built in `main.py`:

```python
from spinh.dft import CollinearVaspJobs, CollinearOrcaJobs, SiteJobs, SaxisJobs, recipes

jobs = CollinearVaspJobs(configs=["uuu", "udu", ...], magmoms=[...], effective_spins=[...])
jobs = CollinearOrcaJobs(configs=[...], spins=[...], effective_spins=[...], expected_moments=[...],
                         n_cpus=32, mem_per_core=16000)
jobs = SiteJobs(magnetic_moments=[7, 7, 4], directions_of="magmom")    # or "spin"
jobs = SaxisJobs(wave_dir="0", e_ref=-1834.36, max_energy=np.inf, de0=1e-8)
```

Sites are numbered from 0 in `main.py` and from 1 in the fit inputs (`fit.yaml`).

## Directions and frames (`site_vasp`, `saxis_vasp`)

Directions are (θ, φ) in degrees. A direction is named by its index on a grid of 1801 θ × 3600 φ
values (0.1°), both poles being single points, so a configuration's directory name is fixed by
its directions.

Configurations are stored in the global frame. Each site can have a local frame:

```python
jobs.set_frames({0: (42.9, 102.5), 1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]]})   # SiteJobs
jobs.set_frame((98.3, 330.1))                                                # SaxisJobs
```

A frame is given by its z direction (θ, φ), with the global x axis in its xz plane, or by a 3×3
matrix whose rows are its axes. Sites not listed use the global frame. `frame="local"` in the
methods below means directions in these frames. The local directions printed and written to the
tables are always computed from the current frames, so changing the frames never leaves stale
local directions behind.

| Method | |
|---|---|
| `add([[(θ1, φ1), (θ2, φ2), …], …], frame="global")` | explicit configurations, one direction per site |
| `add(recipes.updown(["uud", …]))` | collinear configurations (u: θ = 0, d: θ = 180) |
| `add_uniform(directions, frame)` | for each direction, every site along it |
| `add_exchange_pair(i, j, rest=(180, 0), frame="local")` | the 18 configurations for the exchange matrix of sites i, j (`SiteJobs`; see below) |
| `set_direction(site, (θ, φ), frame)` | the same direction of one site in every configuration |
| `flip(site)` | reverse one site in every configuration |

Directions for `add_uniform`: `recipes.meridian(phi, theta_min, theta_max, n)`,
`recipes.parallel(theta, phi_min, phi_max, n)`, and `recipes.great_circles(n_theta=7, n_phi=12)`
(the meridians φ = 0, 180, 90, 270 and the equator). A configuration whose directory exists
already is skipped.

`SiteJobs` takes the directions of the magnetic moments (`directions_of="magmom"`), or of the
spins (`directions_of="spin"`, moments opposite).

### Exchange matrices

`jobs.add_exchange_pair(i, j)` adds, for every axis a of site i and b of site j (default x, y, z
of their frames), the configuration with site i along −a and the one with site i along +a, site j
along b and the other sites along `rest`. Each such pair of configurations is one equation of
`spinh fit anisotropic-exchange`. `python main.py energies` then also writes
`exchange_energies.yaml`, which `fit.yaml` reads with `energies: exchange_energies.yaml`.

## Actions

| Action | |
|---|---|
| `print`, `print-global` | configurations (local frames / global frame) and their directories |
| `setup` | create the directories, write the inputs, link POTCAR (`saxis_vasp`: WAVECAR from the job `wave_dir`) |
| `submit` | `sbatch` in every directory |
| `check`, `resubmit` | convergence; resubmit the jobs that are not converged |
| `energies` | `energies.dat` (`site_vasp` with exchange pairs: also `exchange_energies.yaml`) |
| `all-energies` | add every job directory found, then `energies` |
| `restart-test`, `restart` | resubmit unfinished jobs, optionally with WAVECAR/CHGCAR from the finished job with the nearest directions (`jobs.restart_options`) |
| `magmoms` | site moments (`collinear_vasp`), their angles to the initial directions (`site_vasp`), or the Loewdin populations checked against `expected_moments` (`collinear_orca`, `[tolerance]`) |
| `dipoles`, `dipoles-crystal` | `dipoles.dat` (dipole correction, or Berry-phase dipole), `site_vasp` |
| `site-moments <i>`, `occmat` | spin and orbital moments of site i, occupation matrices, `site_vasp` |
| `backup`, `restart-inputs` | `collinear_orca`: archive unconverged runs as `run_N.*`; inputs that restart from `dft.gbw` |

Convergence: `site_vasp` and `collinear_vasp` look for "EDIFF is reached" in OUTCAR;
`saxis_vasp` takes the energy change in the job's `output` (below `de0`) and keeps only
energies within `max_energy` of `e_ref`; `collinear_orca` counts "SCF CONVERGED" in `orca.out`.

## Output files

`energies.dat` and `dipoles.dat` of the noncollinear jobs have one line per configuration: θ, φ
of every site, then the energy (eV) or px, py, pz, |p| (e Å). The block in the global frame comes
first and is commented out, so `np.loadtxt` reads the block in the local frames.

`energies.dat` of the collinear jobs has the signed effective spin of each site, then the energy;
its first line names the unit (eV for VASP, Ha for ORCA). Jobs without an energy are left out.

## Job directories from before the plain input files

Their inputs were Python strings in `data.py`. `spinh migrate <dir>` writes them as `INCAR`,
`POSCAR`, `KPOINTS`, `vasp.job` (or `orca_hs.inp`, `orca_bs.inp`, `orca.job`), never overwriting
a file, and prints the other settings of `data.py` with the argument they now go to. Then
rewrite `main.py` as in the template.
