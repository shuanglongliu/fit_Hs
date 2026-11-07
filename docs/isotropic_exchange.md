# Isotropic exchange

$$H = \text{factor} \sum_k J_k \sum_{\langle ij \rangle \in k} \mathbf{S}_i \cdot \mathbf{S}_j$$

`factor` is −2 by default.

DFT: `collinear_vasp` or `collinear_orca` (`python main.py energies` writes `energies.dat`).
Fit: `spinh new fit_isotropic_exchange <dir>`, then in it

```bash
spinh fit isotropic-exchange
```

`fit.yaml`:

```yaml
data: energies.dat     # signed effective spins of every site, then the energy;
                       # the first line names the unit (eV or Ha)
factor: -2
threshold: null        # cm^-1, optional
pairs:                 # J label: the 1-indexed site pairs sharing that J
  J1: [[1, 2], [3, 4]]
  J2: [[1, 3]]
```

Energy differences are taken relative to the last configuration. With as many differences as
*J* values the system is solved exactly, with more by least squares (the residual is printed);
fewer is an error. `threshold` drops configurations more than that many cm⁻¹ above the lowest.

Outputs: the *J* values in cm⁻¹ and meV on screen, `Js.csv` (cm⁻¹) and `fitting_quality.pdf`.
