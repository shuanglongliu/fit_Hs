# Anisotropic exchange

$$H = \text{factor} \sum_{\langle ij \rangle} \mathbf{S}_i^{\mathsf T}\, \mathbf{J}_{ij}\, \mathbf{S}_j$$

Each **J**<sub>*ij*</sub> is written in a reference frame attached to the pair. Every energy
pair (conf1, conf2) gives one equation

$$E_2 - E_1 = \text{factor} \sum_{\langle ij \rangle} \sum_{\alpha\beta} J_{ij}^{\alpha\beta}
\left( S_{i\alpha}^{(2)} S_{j\beta}^{(2)} - S_{i\alpha}^{(1)} S_{j\beta}^{(1)} \right)$$

with the spins rotated into the frame of pair *ij*. With *n* pairs there are 9*n* unknowns,
ordered pair-1 xx, xy, …, pair-n zz.

DFT: `site_vasp`, with `jobs.add_exchange_pair(i, j)` in its `main.py` for every pair (see
[dft_jobs.md](dft_jobs.md)); `python main.py energies` then writes `exchange_energies.yaml`.
Fit: `spinh new fit_anisotropic_exchange <dir>`, then in it

```bash
spinh fit anisotropic-exchange
```

*J* is printed in meV and cm⁻¹, per pair and then on one line per unit.
`quantum_correction: true` multiplies by 2*S*<sub>i</sub>*S*<sub>j</sub> / (2*S*<sub>i</sub>*S*<sub>j</sub> + min(*S*<sub>i</sub>, *S*<sub>j</sub>)) (Ruiz et al. 2003).

## `fit.yaml`

```yaml
factor: -2
quantum_correction: false
spins: [1.5, 1.5, 1.5]

exchange:              # one entry per pair: the unknown J matrices
  - pair: [1, 2]
    reference_frame: [ex_x, ex_y, ex_z,  ey_x, ey_y, ey_z,  ez_x, ez_y, ez_z]

energies:              # one entry per equation
  - pair: [1, 2]       # frame in which the angles are given
    conf1: [theta_1, phi_1, theta_2, phi_2, ..., energy_eV]
    conf2: [theta_1, phi_1, theta_2, phi_2, ..., energy_eV]
```

or `energies: exchange_energies.yaml`, the file written by `python main.py energies`. An entry with
`frame: global` has its angles in the global frame (as in that file); otherwise they are in the
frame of its `pair`.

At least 9 energy pairs per exchange pair are needed; with more the system is solved by least
squares and the residual is printed.

## Helpers

| Command | |
|---|---|
| `spinh tool pair-frame I J --structure POSCAR` | frame of pair (I, J): e<sub>x</sub> along the bond, e<sub>z</sub> the part of the global z axis normal to it (needs ase); prints a `reference_frame:` line |
| `spinh tool average-halves` | `energies_average.dat`: energy of row k averaged with row k + m/2 of `energies.dat` |
| `spinh tool entries I J [--n-pair 9]` | `energies:` entries for sites I, J from `energies.dat`, row k paired with row k + 9 (for job sets built without `add_exchange_pair`) |

To express a fitted **J** of a pair in the global frame, or that of a pair related by a rotation
*R*: `spinh.models.anisotropic_exchange.J_local_to_global(emat, J, R)` (the pair frame is
(*R* emat<sup>T</sup>)<sup>T</sup>), e.g. for C₃-related pairs
`R = spinh.geometry.euler_zyz(120, 0, 0)`.
