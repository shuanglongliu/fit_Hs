"""Full 3x3 exchange matrices from noncollinear DFT energies.

    H = factor * sum_{ij} S_i^T J_ij S_j

Each J_ij is written in a reference frame attached to the pair. Every energy pair (conf1, conf2)
in fit.yaml gives one equation

    E2 - E1 = factor * sum_{ij} sum_{ab} J_ij^{ab} (S2_ia S2_jb - S1_ia S1_jb)

with the spin components rotated into the frame of pair ij. The angles of an energy pair are given
in the frame of its pair, or in the global frame with "frame: global" (as written by the
"energies" action of SiteJobs.add_exchange_pair). Energies are read in eV; J is reported in meV
and cm^-1.
"""

import numpy as np
import yaml

from spinh import linfit, units
from spinh.geometry import sph2cart

LABELS = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]


def read_input(cfg):
    """spins, spin_pairs, reference_frames (columns = local axes) and energy_pairs from the
    settings; ``energies`` is a list of entries or the name of a YAML file holding one."""
    energies = cfg["energies"]
    if isinstance(energies, str):
        with open(energies) as f:
            energies = yaml.safe_load(f)
    return {
        "spins": np.array(cfg["spins"], dtype=float),
        "spin_pairs": [list(entry["pair"]) for entry in cfg["exchange"]],
        # reference_frame lists ex, ey, ez one after another; transpose to make them columns
        "reference_frames": [np.array(entry["reference_frame"], dtype=float).reshape(3, 3).T
                             for entry in cfg["exchange"]],
        "energy_pairs": energies,
    }


def spin_vectors(conf, spins):
    """Spin vectors (n_site, 3) of one configuration [theta_1, phi_1, ..., energy]."""
    return np.array([sph2cart([S, conf[2 * i], conf[2 * i + 1]]) for i, S in enumerate(spins)])


def coefficient_matrix(inp, factor):
    """C (n_energy_pair, 9 n_spin_pair) and dE = E2 - E1 in meV."""
    spins, spin_pairs, frames = inp["spins"], inp["spin_pairs"], inp["reference_frames"]
    n_site = len(spins)
    C = np.zeros((len(inp["energy_pairs"]), 9 * len(spin_pairs)))
    dE = np.zeros(len(inp["energy_pairs"]))
    for row, energy_pair in enumerate(inp["energy_pairs"]):
        # The angles are given in the frame of the listed pair, or with "frame: global" in the global frame
        if energy_pair.get("frame", "pair") == "global":
            frame_from = np.eye(3)
        else:
            frame_from = frames[spin_pairs.index(list(energy_pair["pair"]))]
        S1 = spin_vectors(energy_pair["conf1"], spins)
        S2 = spin_vectors(energy_pair["conf2"], spins)
        for k, (i, j) in enumerate(spin_pairs):
            T = frames[k].T @ frame_from        # into the frame of pair k
            S1_k, S2_k = S1 @ T.T, S2 @ T.T
            C[row, 9 * k:9 * (k + 1)] = factor * (np.outer(S2_k[i - 1], S2_k[j - 1])
                                                  - np.outer(S1_k[i - 1], S1_k[j - 1])).flatten()
        dE[row] = units.MEV_PER_EV * (energy_pair["conf2"][2 * n_site] - energy_pair["conf1"][2 * n_site])
    return C, dE


def get_Js(inp, factor=-2.0, quantum_correction=False):
    """J components in meV, ordered pair-1 xx, xy, ..., pair-n zz."""
    C, dE = coefficient_matrix(inp, factor)
    Js = linfit.solve(C, dE, require_determined=True)
    linfit.report(C, Js, dE, "meV")
    if quantum_correction:
        # Empirical correction of Ruiz et al. (2003)
        for k, (i, j) in enumerate(inp["spin_pairs"]):
            Si, Sj = inp["spins"][i - 1], inp["spins"][j - 1]
            Js[9 * k:9 * (k + 1)] *= 2 * Si * Sj / (2 * Si * Sj + min(Si, Sj))
    return Js


def print_Js(Js, spin_pairs):
    """One block per pair (meV and cm^-1), then all components on one line per unit."""
    for k, (i, j) in enumerate(spin_pairs):
        print(f"Pair {i} {j}{'meV':>15s}{'cm^-1':>10s}")
        for label, J in zip(LABELS, Js[9 * k:9 * (k + 1)]):
            print(f"  J_{label} = {J:10.3f}{J * units.CM_PER_MEV:10.3f}")
    print("J (meV):  " + "".join(f" {J:9.3f}" for J in Js))
    print("J (cm^-1):" + "".join(f" {J * units.CM_PER_MEV:9.3f}" for J in Js))


def fit(cfg):
    """Settings (fit.yaml): spins, exchange, energies, factor (default -2), quantum_correction."""
    inp = read_input(cfg)
    Js = get_Js(inp, float(cfg.get("factor", -2.0)), bool(cfg.get("quantum_correction", False)))
    print_Js(Js, inp["spin_pairs"])
    return Js


def average_halves(fin="energies.dat", fout="energies_average.dat"):
    """Average the energy (last column) of row i and row i + n/2: a configuration and its
    counterpart listed in the second half of the file (e.g. with the other spins reversed)."""
    data = np.loadtxt(fin)
    m, n = data.shape
    half = m // 2
    data[:half, n - 1] = (data[:half, n - 1] + data[half:m, n - 1]) / 2
    np.savetxt(fout, data[:half], fmt="%15.8f")


def energy_pair_entries(i, j, fin="energies.dat", n_pair=9):
    """Print fit.yaml "energies" entries for sites i and j: row k is paired with row k + n_pair."""
    data = np.loadtxt(fin)
    ii, jj = min(i, j), max(i, j)
    cols = [2 * (ii - 1), 2 * ii - 1, 2 * (jj - 1), 2 * jj - 1]
    for k in range(n_pair):
        conf1 = list(data[k, cols]) + [data[k, -1]]
        conf2 = list(data[k + n_pair, cols]) + [data[k + n_pair, -1]]
        print(f"\n  - pair: [{i}, {j}]")
        print("    conf1:   [{:6.1f}, {:6.1f}, {:6.1f}, {:6.1f}, {:18.9f}]".format(*conf1))
        print("    conf2:   [{:6.1f}, {:6.1f}, {:6.1f}, {:6.1f}, {:18.9f}]".format(*conf2))


def pair_frame(structure, i, j):
    """Frame of pair (i, j): ex along the bond, ez the part of the global z axis normal to it.

    Prints it as a fit.yaml reference_frame entry.
    """
    from ase.io import read

    positions = read(structure).get_positions()
    ex = positions[j - 1] - positions[i - 1]
    ex /= np.linalg.norm(ex)
    ez = np.array([0., 0., 1.]) - ex[2] * ex
    if np.linalg.norm(ez) < 1.e-8:
        raise ValueError(f"Pair {i} {j}: bond is parallel to the global z axis.")
    ez /= np.linalg.norm(ez)
    emat = np.array([ex, np.cross(ez, ex), ez])
    print(f"Pair {i:3d} {j:3d}")
    print("reference_frame: [" + ", ".join(f"{v:15.8f}" for v in emat.flatten()) + "]")
    return emat


def J_local_to_global(emat, J, rotmat=np.eye(3)):
    """J of a pair whose frame is (rotmat emat^T)^T, expressed in the global frame."""
    emat_ = (rotmat @ np.transpose(emat)).T
    return emat_.T @ J @ emat_
