# Single-ion anisotropy

DFT: `saxis_vasp` (the whole system's spin axis) or `site_vasp` (one site rotated, the others
fixed). `python main.py energies` writes `energies.dat`.

Analysis: `spinh new fit_anisotropy <dir>` gives `anisotropy.yaml` (one section per task) and
synthetic example data. In the directory:

| `spinh anisotropy …` | |
|---|---|
| `average` | average E(**n**) and E(−**n**) (removes the exchange with fixed neighbours) → `<input>_average.dat` |
| `change-frame` | directions in a frame with z along `z_direction`, or the lowest/highest-energy direction → `new_frame.dat` |
| `mae` | out-of-plane and in-plane anisotropy energies from θ = 0 and the θ = 90° circle |
| `plot` | polar map of the energies → `<input>.pdf` |
| `fit-DE` | *D* [*S*<sub>z</sub>² − *S*(*S* + 1)/3] + *E* (*S*<sub>x</sub>² − *S*<sub>y</sub>²) + *E*<sub>0</sub> in a fitted eigenframe; quantum or classical spins → `energies_fitted.dat`, `d_spin-orbit.dat`, `fitting_quality.pdf` |
| `fit-D-tensor` | **S**<sup>T</sup> **D** **S** + *E*<sub>0</sub> with a general symmetric **D**, then *D*, *E* and the eigenframe |
| `fit-curve` | *a* (θ − *b*)² + *c* along one cut |
| `fit-surface` | *a* θ′² + *c*, θ′ measured from a fitted axis |

`input` names the table (`<input>.dat`) and `columns` the 1-indexed θ, φ and energy columns.
Energies are read in eV and fitted in meV relative to the minimum; zero energies (unfinished jobs)
are ignored. `spin` is the spin quantum number of the quantum models.

## Fit parameters

`bounds` gives `[min, max]` for each parameter, or a single number to fix it; `p0` gives start
values (default: the middle of the bounds). Angles are in degrees.

The eigenframe is the frame of the data (`frame`, default the identity) rotated by the extrinsic
zyz Euler angles (α, β, γ) — `spinh.geometry.rotate_frame` — the same convention as the Euler
angles of the polarization constraints:

- the eigenframe's z axis points along (θ, φ) = (β, γ);
- α turns its x and y axes about that z axis.

So `fit-DE` with β and γ free finds an easy axis in any direction and α orients *E*. Fixing
γ = 0 and β ≈ 0 (the template) fits only the in-plane orientation α about a known z axis.
`fit-surface` fits only the axis (β, γ).
