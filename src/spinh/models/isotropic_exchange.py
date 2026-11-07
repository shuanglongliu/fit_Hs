"""Isotropic exchange constants from collinear DFT energies.

    H = factor * sum_k J_k sum_{ij in k} S_i . S_j

Inputs: energies.dat (signed effective spins of each site, then the energy; the first line names
the unit, eV or Ha) and pairs.yaml (J label -> list of 1-indexed site pairs).
Energy differences are taken relative to the last configuration. Output: Js.csv (cm^-1) and
fitting_quality.pdf.
"""

import numpy as np
import yaml

from spinh import linfit, units
from spinh.plotting import plot_fit


def read_energies(fname="energies.dat", energy_threshold=None):
    """Spins (n_conf, n_site) and energies relative to the first configuration in cm^-1.

    Configurations more than ``energy_threshold`` (cm^-1) above the lowest one are dropped.
    """
    with open(fname) as f:
        header = f.readline()
    if "eV" in header:
        unit_factor = units.CM_PER_EV
    elif "Ha" in header:
        unit_factor = units.CM_PER_HA
    else:
        raise ValueError(f"Unknown energy unit in {fname} (expected 'eV' or 'Ha').")

    data = np.loadtxt(fname, ndmin=2)
    data[:, -1] = (data[:, -1] - data[0, -1]) * unit_factor
    if energy_threshold is not None:
        mask = data[:, -1] - data[:, -1].min() <= energy_threshold
        if not mask.all():
            print(f"Excluding {np.sum(~mask)} configuration(s) with energy > {energy_threshold} cm^-1 "
                  f"above the minimum.")
            data = data[mask]
    return data[:, :-1], data[:, -1]


def get_sumSS(spins, pairs_dict):
    """sumSS[k, c] = sum over the pairs of J_k of S_i S_j in configuration c."""
    return np.array([sum(spins[:, i - 1] * spins[:, j - 1] for i, j in pairs)
                     for pairs in pairs_dict.values()])


def fit(dirname=".", factor=-2.0, energy_threshold=None):
    """Fit the J values (cm^-1); writes Js.csv and fitting_quality.pdf in the working directory."""
    spins, energies = read_energies("energies.dat", energy_threshold)
    with open("pairs.yaml") as f:
        pairs_dict = yaml.safe_load(f)

    sumSS = get_sumSS(spins, pairs_dict)
    C = factor * (sumSS[:, :-1] - sumSS[:, [-1]]).T
    dE = energies[:-1] - energies[-1]
    Js = linfit.solve(C, dE, require_determined=True)
    linfit.report(C, Js, dE, "cm^-1")

    with open("Js.csv", "w") as f:
        f.write(f"# Js in cm^-1, factor = {factor:g}\n")
        f.write(",".join(pairs_dict) + "\n")
        f.write(",".join(f"{J:.2f}" for J in Js) + "\n")
    print(f"{'':6s} {'cm^-1':>10s} {'meV':>10s}")
    for label, J in zip(pairs_dict, Js):
        print(f"{label:6s} {J:10.2f} {J / units.CM_PER_MEV:10.3f}")
    print("Saved Js.csv (cm^-1).")
    plot_fit(dE, C @ Js, "fitting_quality.pdf", "configuration", r"$\Delta E$ (cm$^{-1}$)")
    return Js
