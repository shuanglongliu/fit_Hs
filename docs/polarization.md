# Spin-driven electric polarization

DFT: `site_vasp` (`python main.py dipoles` writes `dipoles.dat`). Fit: `spinh new
fit_polarization_<model> <dir>`, then in it `spinh polarization <model>`.

All models fit *y*, the change of the dipole of every configuration relative to the first one,
by minimizing |*A x* − *y*|². When the data do not determine every parameter, the minimum-norm
solution is returned with a warning that says how many combinations are undetermined.

`dipoles.dat`: θ, φ (degrees) of every site, then *p*<sub>x</sub>, *p*<sub>y</sub>,
*p*<sub>z</sub> (e Å); the first row is the reference.

## `alpha-pi`

$$\mathbf{p} = \sum_{\langle ij \rangle} \Big[ \alpha_{ij} \big( \mathbf{S}_i (\mathbf{S}_j \cdot \mathbf{e}_{ij}) - (\mathbf{S}_i \cdot \mathbf{e}_{ij}) \mathbf{S}_j \big) + \boldsymbol{\pi}_{ij} (\mathbf{S}_i \cdot \mathbf{S}_j) \Big]$$

**e**<sub>ij</sub> is the unit vector from site *i* to site *j*. `input.yaml`: `spins`, `pairs`
(1-indexed), `positions` (Å), `constraints_alpha` (`[pair_i, pair_j]`: equal α) and
`constraints_pi` (`[pair_i, pair_j, α, β, γ]`: **π**(pair_j) = *R* **π**(pair_i), *R* from
extrinsic zyz Euler angles in degrees). Output: `log`, `fitting_quality.pdf`.

## `full-p-v1` and `full-p-v2`

$$p_\gamma = \sum_{\langle ij \rangle} \sum_{\alpha\beta} S_{i\alpha}\, P^{ij}_{\gamma\alpha\beta}\, S_{j\beta}$$

stored as `x[27 k + 9 γ + 3 α + β]`. The two versions solve the same problem and give the same
result:

- `full-p-v1`: every pair is unknown; each `constraints` entry `[pair_i, pair_j, α, β, γ]`
  imposes *P*<sup>pair_j</sup> = (*R* ⊗ *R* ⊗ *R*) *P*<sup>pair_i</sup>. Output: `P_total.log`,
  `fit.log`, `fitting_quality.pdf`.
- `full-p-v2`: only the pairs in `irr_pairs` are unknown; each `transformations` entry builds
  pair_j from the irreducible pair_i. `emats` gives the local frame of each irreducible pair
  (e<sub>x</sub>, e<sub>y</sub>, e<sub>z</sub>). Output: `fit.log` (with percentage deviations),
  `fitting_quality.pdf`, `P_total`, `P_J`, `P_D`, `P_Gamma` and `pi_and_d` logs in the global frame
  and in the local frames (`*_loc.log`), `relative_contribution.log`, and `P_for_latex.log`
  (local frames, 10⁻⁵ e Å).

*P*<sup>J</sup>, *P*<sup>D</sup>, *P*<sup>Γ</sup> are the isotropic, antisymmetric and anisotropic
symmetric parts over (α, β); **π**<sub>γ</sub> = *P*<sup>J</sup><sub>γxx</sub> and
**d**<sub>γ</sub> = (*P*<sup>D</sup><sub>γyz</sub>, *P*<sup>D</sup><sub>γzx</sub>, *P*<sup>D</sup><sub>γxy</sub>).

With the example data the problem is underdetermined (6 of 54 combinations), so more spin
configurations are needed to determine every element. The earlier Julia version 1 solved a
singular KKT system that added an arbitrary null-space component, which is why the two versions
used to disagree while giving the same fitting errors.
