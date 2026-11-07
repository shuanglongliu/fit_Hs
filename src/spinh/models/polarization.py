"""Spin-driven electric polarization of a magnetic cluster, fitted to DFT dipoles.

Three models, all fitted to y = the change of the dipole of every configuration relative to the
first one (input.yaml and dipoles.dat in the working directory; dipoles.dat has theta, phi of
every site in degrees, then px, py, pz in e Ang):

alpha-pi    p = sum_ij [alpha_ij (S_i (S_j . e_ij) - (S_i . e_ij) S_j) + pi_ij (S_i . S_j)]
            x = [alpha_1..alpha_n, pi_1x, pi_1y, pi_1z, ..., pi_nz]
full-p-v1   p_g = sum_ij sum_ab S_ia P^ij_gab S_jb, every pair unknown, symmetry as constraints
full-p-v2   the same model with only the irreducible pairs unknown; the other pairs by rotation

Full tensors are stored as x[27 k + 9 gamma + 3 alpha + beta].
"""

import numpy as np
import yaml

from spinh import linfit
from spinh.geometry import euler_zyz, sph2cart
from spinh.plotting import plot_fit

UNIT = "e Ang"


def read_input():
    with open("input.yaml") as f:
        return yaml.safe_load(f)


def read_dipoles(spins):
    """Spin vectors (n_conf, n_site, 3) and dipoles (n_conf, 3) from dipoles.dat."""
    data = np.loadtxt("dipoles.dat", ndmin=2)
    n_site = len(spins)
    S = np.array([[sph2cart([spins[j], row[2 * j], row[2 * j + 1]]) for j in range(n_site)] for row in data])
    return S, data[:, 2 * n_site:2 * n_site + 3]


def differences(B, i_conf_ref=0):
    """Rows of B minus the reference row, without the reference, flattened over (conf, gamma)."""
    return np.delete(B - B[i_conf_ref], i_conf_ref, axis=0).reshape(-1, *B.shape[2:])


def deviations(x, y, A, fname, mode="w", percentages=False):
    """Write the deviations of the fit and plot DFT vs fit."""
    y_fit = A @ x
    dy = np.abs(y_fit - y)
    with open(fname, mode) as f:
        if mode == "a":
            f.write("\n\n")
        f.write(f"Maximal deviation: {dy.max():12.6f} {UNIT}.\n")
        f.write(f"Average deviation: {dy.mean():12.6f} {UNIT}.\n")
        if percentages:
            # The first element is always included, the rest only where y >= 0.001
            keep = np.concatenate([[True], y[1:] >= 0.001])
            with np.errstate(divide="ignore", invalid="ignore"):
                dy_per = 100 * dy[keep] / y[keep]
            f.write(f" std of deviation: {dy.std(ddof=1):12.6f} {UNIT}.\n\n")
            f.write(f"Maximal deviation in percentage: {dy_per.max():12.6f} \n")
            f.write(f"Average deviation in percentage: {dy_per.mean():12.6f} \n")
            f.write(f" std of deviation in percentage: {dy_per.std(ddof=1):12.6f} \n")
    plot_fit(y, y_fit, "fitting_quality.pdf", r"Index of $\delta P_\alpha$ for all spin configurations",
             r"$\delta P_\alpha$ (e $\AA$)")


def write_P(f, x, k):
    """The 3x3 matrices P_gamma of block k of x."""
    for gamma, direction in enumerate("xyz"):
        f.write(f"The {direction} component of coupling matrix\n\n")
        for alpha in range(3):
            index = 27 * k + 9 * gamma + 3 * alpha
            f.write("{:12.6f} {:12.6f} {:12.6f}\n".format(*x[index:index + 3]))
        f.write("\n")


def save_P(x, pairs, fname):
    """P of every listed pair (block k of x belongs to pairs[k])."""
    with open(fname, "w") as f:
        for k, (i, j) in enumerate(pairs):
            f.write(f"Pair {i:3d} {j:3d}\n\n")
            write_P(f, x, k)
            f.write("\n")


# ---------------------------------------------------------------------------------------------
# alpha-pi
# ---------------------------------------------------------------------------------------------

def fit_alpha_pi():
    """Writes log and fitting_quality.pdf."""
    inp = read_input()
    spins = np.array(inp["spins"], dtype=float)
    pairs = inp["pairs"]
    positions = np.array(inp["positions"], dtype=float)
    eij = np.array([positions[j - 1] - positions[i - 1] for i, j in pairs])
    eij /= np.linalg.norm(eij, axis=1)[:, None]
    S, P = read_dipoles(spins)
    n_conf, n_pair = S.shape[0], len(pairs)

    B = np.zeros((n_conf, 3, 4 * n_pair))
    for k, (i, j) in enumerate(pairs):
        Si, Sj = S[:, i - 1], S[:, j - 1]
        B[:, :, k] = Si * (Sj @ eij[k])[:, None] - (Si @ eij[k])[:, None] * Sj
        for a in range(3):
            B[:, a, n_pair + 3 * k + a] = np.einsum("ca,ca->c", Si, Sj)
    A, y = differences(B).reshape(-1, 4 * n_pair), differences(P).flatten()

    # alpha(pair_i) = alpha(pair_j); pi(pair_j) = R pi(pair_i)
    constraints_alpha = inp.get("constraints_alpha") or []
    constraints_pi = inp.get("constraints_pi") or []
    C = np.zeros((len(constraints_alpha) + 3 * len(constraints_pi), 4 * n_pair))
    for row, (i_pair, j_pair) in enumerate(constraints_alpha):
        C[row, i_pair - 1], C[row, j_pair - 1] = 1.0, -1.0
    row = len(constraints_alpha)
    for i_pair, j_pair, alpha, beta, gamma in constraints_pi:
        rot = euler_zyz(alpha, beta, gamma)
        for a in range(3):
            C[row, n_pair + 3 * (i_pair - 1):n_pair + 3 * i_pair] = rot[a]
            C[row, n_pair + 3 * (j_pair - 1) + a] = -1.0
            row += 1
    x = linfit.solve(A, y, C, np.zeros(len(C)))

    with open("log", "w") as f:
        for k, (i, j) in enumerate(pairs):
            pi = x[n_pair + 3 * k:n_pair + 3 * k + 3]
            f.write(f"Pair {i:3d} {j:3d}\n\n")
            f.write(f"alpha = {x[k]:12.6f} {UNIT}\n")
            f.write(f"pi = ({pi[0]:12.6f} {pi[1]:12.6f} {pi[2]:12.6f}) {UNIT}\n\n")
    deviations(x, y, A, "log", mode="a")
    return x


# ---------------------------------------------------------------------------------------------
# full P tensor
# ---------------------------------------------------------------------------------------------

def tensor_rotation(alpha, beta, gamma):
    """R x R x R acting on the 27 elements of one P tensor."""
    rot = euler_zyz(alpha, beta, gamma)
    return np.kron(rot, np.kron(rot, rot))


def fit_full_p_v1():
    """All pairs unknown; each constraint [pair_i, pair_j, a, b, g] imposes P(pair_j) = R P(pair_i).
    Writes P_total.log, fit.log and fitting_quality.pdf."""
    inp = read_input()
    pairs = inp["pairs"]
    S, P = read_dipoles(np.array(inp["spins"], dtype=float))
    n_conf, n_pair = S.shape[0], len(pairs)

    B = np.zeros((n_conf, 3, 27 * n_pair))
    for k, (i, j) in enumerate(pairs):
        SiSj = np.einsum("ca,cb->cab", S[:, i - 1], S[:, j - 1]).reshape(n_conf, 9)
        for gamma in range(3):
            B[:, gamma, 27 * k + 9 * gamma:27 * k + 9 * (gamma + 1)] = SiSj
    A, y = differences(B).reshape(-1, 27 * n_pair), differences(P).flatten()

    constraints = inp.get("constraints") or []
    C = np.zeros((27 * len(constraints), 27 * n_pair))
    for ic, (i_pair, j_pair, alpha, beta, gamma) in enumerate(constraints):
        rows = slice(27 * ic, 27 * (ic + 1))
        C[rows, 27 * (i_pair - 1):27 * i_pair] = tensor_rotation(alpha, beta, gamma)
        C[rows, 27 * (j_pair - 1):27 * j_pair] = -np.eye(27)
    x = linfit.solve(A, y, C, np.zeros(len(C)))

    save_P(x, pairs, "P_total.log")
    deviations(x, y, A, "fit.log")
    return x


def decompose_P(x):
    """Isotropic (J), antisymmetric (D) and anisotropic symmetric (Gamma) parts over (alpha, beta)."""
    P = x.reshape(-1, 3, 3, 3)
    x_J = np.einsum("kg,ab->kgab", np.trace(P, axis1=2, axis2=3) / 3, np.eye(3)).flatten()
    x_D = ((P - P.transpose(0, 1, 3, 2)) / 2).flatten()
    return x_J, x_D, x - x_J - x_D


def save_pi_and_d(x_J, x_D, pairs, fname):
    """pi_gamma = P^J_{gamma xx} and d_gamma = (P^D_{gamma yz}, P^D_{gamma zx}, P^D_{gamma xy})."""
    P_J, P_D = x_J.reshape(-1, 3, 3, 3), x_D.reshape(-1, 3, 3, 3)
    with open(fname, "w") as f:
        for k, (i, j) in enumerate(pairs):
            f.write(f"Pair {i:3d} {j:3d}\n\n")
            f.write("pi = {:12.6f} {:12.6f} {:12.6f}\n\n".format(*P_J[k, :, 0, 0]))
            for gamma, direction in enumerate("xyz"):
                d = P_D[k, gamma]
                f.write(f"d{direction} = {d[1, 2]:12.6f} {d[2, 0]:12.6f} {d[0, 1]:12.6f}\n")
            f.write("\n\n")


def save_results(x, pairs, modifier):
    save_P(x, pairs, "P_total" + modifier + ".log")
    x_J, x_D, x_Gamma = decompose_P(x)
    save_P(x_J, pairs, "P_J" + modifier + ".log")
    save_P(x_D, pairs, "P_D" + modifier + ".log")
    save_P(x_Gamma, pairs, "P_Gamma" + modifier + ".log")
    save_pi_and_d(x_J, x_D, pairs, "pi_and_d" + modifier + ".log")


def to_local_frames(x, emats):
    """x_loc_gab = sum emat_gg' emat_aa' emat_bb' x_g'a'b' (rows of emat = local axes)."""
    P = x.reshape(-1, 3, 3, 3)
    return np.array([np.einsum("gh,ai,bj,hij->gab", e, e, e, P[k]) for k, e in enumerate(emats)]).flatten()


def save_relative_contribution(x, A, pairs):
    """Share of each irreducible pair in |dP|, and of its J, D and Gamma parts."""
    blocks = [slice(27 * k, 27 * (k + 1)) for k in range(len(pairs))]
    y_irr = [np.abs(A[:, b] @ x[b]) for b in blocks]
    y_sum_irr = np.sum(y_irr, axis=0)
    x_J, x_D, x_Gamma = decompose_P(x)

    def stats(f, y_per, tag):
        f.write(f"Maximal contribution from {tag:>8s}: {np.nanmax(y_per):6.1f} percent.\n")
        f.write(f"Average contribution from {tag:>8s}: {np.nanmean(y_per):6.1f} percent.\n")
        f.write(f" std of contribution from {tag:>8s}: {np.nanstd(y_per, ddof=1):6.1f} percent.\n\n")

    with np.errstate(divide="ignore", invalid="ignore"), open("relative_contribution.log", "w") as f:
        for k, b in enumerate(blocks):
            i, j = pairs[k]
            f.write(f"Pair {i:3d} {j:3d}: relative contribution = "
                    f"{np.nanmean(100 * y_irr[k] / y_sum_irr):6.1f} percent.\n\n")
            parts = [np.abs(A[:, b] @ part[b]) for part in (x_J, x_D, x_Gamma)]
            for part, tag in zip(parts, ("P_J", "P_D", "P_Gamma")):
                stats(f, 100 * part / sum(parts), tag)


def save_for_latex(x, pairs):
    """P and its J, D, Gamma parts in units of 1e-5 e Ang, as LaTeX table rows."""
    parts = [("Total ", x)] + list(zip(["The isotropic term", "The antisymmetric term",
                                        "The anisotropic symmetric term"], decompose_P(x)))
    with open("P_for_latex.log", "w") as f:
        for k, (i, j) in enumerate(pairs):
            f.write(f"Pair {i:3d} {j:3d}\n\n")
            for title, part in parts:
                f.write(f"{title}\n\n")
                part = 1e5 * part
                for gamma, direction in enumerate("xyz"):
                    f.write(f"The {direction} component of coupling matrix\n\n")
                    for alpha in range(3):
                        index = 27 * k + 9 * gamma + 3 * alpha
                        end = " \\\\\n" if alpha != 2 else "\n"
                        f.write("{:6.1f} & {:6.1f} &  {:6.1f}".format(*part[index:index + 3]) + end)
                    f.write("\n")
                f.write("\n")


def fit_full_p_v2():
    """Only the pairs in irr_pairs are unknown; a transformation [pair_i, pair_j, a, b, g] builds
    P(pair_j) = R P(pair_i) from the irreducible pair_i. emats gives the local frame of each
    irreducible pair (ex, ey, ez one after another).

    Writes fit.log, fitting_quality.pdf, P_total/P_J/P_D/P_Gamma/pi_and_d logs in the global
    frame and in the local frames (*_loc.log), relative_contribution.log and P_for_latex.log
    (local frames).
    """
    inp = read_input()
    pairs, irr_pairs = inp["pairs"], inp["irr_pairs"]
    emats = [np.array(e, dtype=float).reshape(3, 3) for e in inp["emats"]]
    transformations = {t[1]: t for t in inp.get("transformations") or []}
    S, P = read_dipoles(np.array(inp["spins"], dtype=float))
    n_conf = S.shape[0]

    B = np.zeros((n_conf, 3, 27 * len(irr_pairs)))
    for k_pair, (i_site, j_site) in enumerate(pairs, start=1):
        Si, Sj = S[:, i_site - 1], S[:, j_site - 1]
        if k_pair in irr_pairs:
            rot, i_irrp = np.eye(3), irr_pairs.index(k_pair)
        else:
            i_pair, _, alpha, beta, gamma = transformations[k_pair]
            rot, i_irrp = euler_zyz(alpha, beta, gamma), irr_pairs.index(i_pair)
        # B[c, gamma, (gamma', alpha', beta')] = R_{gamma gamma'} (R^T Si)_{alpha'} (R^T Sj)_{beta'}
        block = np.einsum("gh,ca,cb->cghab", rot, Si @ rot, Sj @ rot).reshape(n_conf, 3, 27)
        B[:, :, 27 * i_irrp:27 * (i_irrp + 1)] += block
    A, y = differences(B).reshape(-1, 27 * len(irr_pairs)), differences(P).flatten()
    x = linfit.solve(A, y)

    irr = [pairs[k - 1] for k in irr_pairs]
    deviations(x, y, A, "fit.log", percentages=True)
    save_results(x, irr, "")
    x_loc = to_local_frames(x, emats)
    save_results(x_loc, irr, "_loc")
    save_relative_contribution(x, A, irr)
    save_for_latex(x_loc, irr)
    return x, x_loc


MODELS = {"alpha-pi": fit_alpha_pi, "full-p-v1": fit_full_p_v1, "full-p-v2": fit_full_p_v2}
