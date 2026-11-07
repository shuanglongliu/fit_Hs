# spinh

Spin Hamiltonian parameters from DFT: set up and harvest VASP/ORCA jobs for chosen spin
configurations, then fit the parameters to the energies or dipoles.

| Parameters | DFT jobs (template) | Fit | Docs |
|---|---|---|---|
| Isotropic exchange *J* | `collinear_vasp`, `collinear_orca` | `spinh fit isotropic-exchange` | [isotropic_exchange.md](docs/isotropic_exchange.md) |
| Exchange matrices **J** (3×3 per pair) | `site_vasp` (`add_exchange_pair`) | `spinh fit anisotropic-exchange` | [anisotropic_exchange.md](docs/anisotropic_exchange.md) |
| Single-ion anisotropy (*D*, *E*, **D** tensor) | `saxis_vasp`, `site_vasp` | `spinh fit anisotropy-DE`, … | [single_ion_anisotropy.md](docs/single_ion_anisotropy.md) |
| Spin-driven electric polarization | `site_vasp` (`dipoles`) | `spinh fit polarization-alpha-pi`, … | [polarization.md](docs/polarization.md) |

The DFT job templates are described in [dft_jobs.md](docs/dft_jobs.md).

## Install

Requirements: Python ≥ 3.8 and pip ≥ 21.3. numpy, scipy, pyyaml and matplotlib are installed
automatically; ase is optional.

```bash
git clone https://github.com/shuanglongliu/spinh.git
cd spinh

# optional: a separate environment
conda create -n spinh python=3.11 && conda activate spinh     # or: python -m venv .venv && source .venv/bin/activate

pip install -e .                  # the spinh package and the spinh command
pip install -e ".[structure]"     # also ase, needed for `spinh tool pair-frame`
```

`-e` (editable) links the installed package to this checkout, so a `git pull` takes effect
without reinstalling; `pip install .` (without `-e`) also works. Job scripts (`main.py`) import
`spinh`, so they run from any directory while the environment with `spinh` is active.

Check the install:

```bash
spinh -h
python -m unittest discover tests     # from the repository root
```

To remove it: `pip uninstall spinh`.

## A new system

```bash
spinh new                         # list the templates
spinh new site_vasp ~/runs/Mn3    # copy one: INCAR, POSCAR, KPOINTS, vasp.job, main.py
cd ~/runs/Mn3
# 1. edit INCAR, POSCAR, KPOINTS, vasp.job; put POTCAR here
# 2. main.py: the spin configurations
python main.py print              # check the configurations and directory names
python main.py setup              # write the inputs
python main.py submit
python main.py check              # later: convergence
python main.py energies           # energies.dat (dipoles: dipoles.dat)
```

then fit, e.g. `spinh new fit_anisotropy ~/runs/Mn3/fit`, point `data:` in its `fit.yaml` at the
energies and run `spinh fit anisotropy-DE`. Every fit reads `fit.yaml`, which names the data
file. `spinh -h`, `spinh <command> -h` and `python main.py help` list the options. Job directories
made before the inputs became plain files: `spinh migrate <dir>` (see
[dft_jobs.md](docs/dft_jobs.md)).

## Layout

```
spinh/
├── geometry.py, sphere.py        # spherical coordinates, frames, Euler rotations; direction grid
├── spin_operators.py, units.py
├── linfit.py, tables.py, plotting.py
├── dft/                          # job sets for `python main.py ACTION`
│   ├── base.py                   #   JobSet: setup, submit, check, resubmit, energies, help
│   ├── collinear.py, orca.py     #   'u'/'d' configurations: VASP, ORCA
│   ├── noncollinear.py           #   per-site directions (SiteJobs), one spin axis (SaxisJobs)
│   ├── recipes.py                #   ready-made directions and configurations
│   └── vasp.py, slurm.py, migrate.py
├── models/                       # isotropic and anisotropic exchange, single-ion anisotropy, polarization
├── templates/                    # starting points for `spinh new`
│   ├── collinear_vasp, collinear_orca, saxis_vasp, site_vasp   # input files + main.py
│   └── fit_*                     # fit.yaml with example data
└── cli.py                        # the spinh command
docs/
tests/                            # python -m unittest discover tests
```

Each subdirectory of the earlier layout (`collinear_exchange/`, `anisotropic_exchange/`,
`magnetic_anisotropy/`, `polarization/`) was its own repository; the full history is in
`git log`.

## License

GNU General Public License v3.0; see [LICENSE](LICENSE).
