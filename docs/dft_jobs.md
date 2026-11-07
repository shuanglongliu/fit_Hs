# DFT jobs

A job directory holds `data.py` and `main.py`, created by `spinh new <template> <dir>`:

- `data.py` is the calculation: input templates (INCAR, POSCAR, KPOINTS or the ORCA input), the
  Slurm job script and system settings. `root_dir` is the directory of `data.py`; the job
  directories are made under it.
- `main.py` lists the spin configurations and runs an action: `python main.py <action>`.
  `python main.py help` lists the actions.

| Template | Spins | INCAR fields | Directory names |
|---|---|---|---|
| `collinear_vasp` | up/down per site (`MAGMOM`, `NUPDOWN`) | `{magmom}`, `{nupdown}`, `{U}`, `{J}` | `uudu…` |
| `collinear_orca` | broken-symmetry ORCA (`FlipSpin`, `FinalMs`) | see `spinh.dft.orca` | `uudu…` |
| `saxis_vasp` | one axis for the whole system (`SAXIS`) | `{sstring}` | `6476401` |
| `site_vasp` | one moment direction per site (`MAGMOM`, `M_CONSTR`, `RWIGS`) | `{mstring}`, `{rwigs}` | `4319101_0_0` |

POTCAR (and, for ORCA, `orca.xyz`) goes in `root_dir`.

## Directions and frames (`saxis_vasp`, `site_vasp`)

Directions are (θ, φ) in degrees. A direction is named by its index on a grid of 1801 θ × 3600 φ
values (0.1°), both poles being single points, so a configuration's directory name is fixed by
its directions. Directions are kept in the global frame and in local frames:

- `site_vasp`: `myjob.setup_local_ref_frames(sites=[0, 1], z_directions=[(θ, φ), …])` gives
  each listed site a frame whose z axis is along (θ, φ) (global x in its xz plane);
  `from_file=True` reads `data.emats_file` instead. `saxis_vasp`:
  `myjob.setup_local_ref_frames(z_direction=(θ, φ))` or `from_file=True` (`data.emat_file`).
- `local_ref_frame=True` in the `add_*` methods means the directions are in the local frames.
- Changing the frames does not change stored configurations: call
  `myjob.transform_configurations(global2local=True)` to recompute their local directions.

`site_vasp` takes the directions of the magnetic moments (`directions_of="magmom"`), or of the
spins (`directions_of="spin"`, moments opposite).

Building configurations (`SiteJobs`; `SaxisJobs` has the one-site equivalents `add_directions`,
`add_thetas`, `add_phis`, `transform_directions`):

| Method | |
|---|---|
| `add_configurations([[[θ1, φ1], [θ2, φ2], …], …], local_ref_frame)` | explicit configurations |
| `add_thetas(phi, theta_min, theta_max, ntheta, local_ref_frame)` | all sites along one meridian (also `add_thetas_colinear_spin`/`_magmom`) |
| `add_phis(theta, phi_min, phi_max, nphi, local_ref_frame)` | all sites along one parallel |
| `set_direction_for_one_site(i_site, direction, local_ref_frame)` | the same direction of one site in every configuration |
| `flip_direction_for_one_site(i_site)` | reverse one site in every configuration |

## Actions

| Action | |
|---|---|
| `print`, `print-global` | configurations in the local / global frames, with their directories |
| `setup` | create the directories, write the inputs, link POTCAR (`saxis_vasp`: WAVECAR from `data.wave_dir`) |
| `submit` | `sbatch vasp.job` in every directory |
| `check`, `resubmit` | convergence; resubmit the unconverged jobs |
| `restart-test`, `restart` | resubmit unfinished jobs, optionally with WAVECAR/CHGCAR from the finished job with the nearest directions (`myjob.restart_options`) |
| `energies`, `all-energies` | `energies.dat`; `all-energies` first adds every job directory found |
| `dipoles`, `dipoles-crystal` | `dipoles.dat` (dipole correction, or Berry-phase dipole) — `site_vasp` |
| `check-magmoms`, `moments <i>`, `occmat` | moment directions, site spin and orbital moments, occupation matrices — `site_vasp` |

Convergence: `site_vasp` and `collinear_vasp` look for "EDIFF is reached" in OUTCAR;
`saxis_vasp` takes the energy change in the job's `output` (below `myjob.de0`) and keeps only
energies within `myjob.max_energy` of `myjob.e_ref`.

`collinear_vasp` actions: `setup`, `submit`, `check`, `resubmit`, `magmoms`, `energies`.
`collinear_orca` actions: `setup`, `submit`, `check`, `backup`, `restart-inputs`, `resubmit`,
`magmoms [tol]`, `energies [tol]`; see the docstring of `spinh.dft.orca`.

## Output files

`energies.dat` and `dipoles.dat` of the noncollinear jobs have one line per configuration: θ, φ
of every site, then the energy (eV) or px, py, pz, |p| (e Å). The block in the global frame comes
first and is commented out, so `np.loadtxt` reads the block in the local frames.

`energies.dat` of the collinear jobs has the signed effective spin of each site, then the energy;
its first line names the unit (eV for VASP, Ha for ORCA).
