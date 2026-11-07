"""Single-ion magnetic anisotropy from the energies of rotated spins.

Every task reads its settings from fit.yaml in the working directory: the top-level keys and the
section named after the task (see templates/fit_anisotropy/fit.yaml). ``data`` is a text table and
``columns`` its 1-indexed (theta, phi, energy) columns, angles in degrees and energies in eV.

Frames: the anisotropy eigenframe is the DFT frame ``frame`` (rows = axes, default identity)
rotated by the extrinsic zyz Euler angles (alpha, beta, gamma), see spinh.geometry.rotate_frame:
its z axis points along (theta, phi) = (beta, gamma) and alpha turns x and y about it.
"""

import os

import numpy as np
from scipy.optimize import curve_fit

from spinh import units
from spinh.geometry import (cart2sph, change_frame_angles, change_frame_sph, get_emat_local,
                          opposite_direction, rotate_frame, sph2cart)
from spinh.plotting import plot_fit
from spinh.spin_operators import expectation, spin_operators, top_state


def load(cfg):
    """The table ``data`` and its 0-indexed columns."""
    data = np.loadtxt(cfg["data"], ndmin=2)
    return data, [c - 1 for c in cfg["columns"]]


def stem(cfg):
    """``data`` without its extension, for the names of derived files."""
    return os.path.splitext(cfg["data"])[0]


def relative_meV(energies, factor=1.0):
    """1000 factor (E - E_min); zeros (unfinished jobs) become NaN."""
    emin = np.nanmin(energies)
    return np.where(energies == 0.0, np.nan, units.MEV_PER_EV * factor * (energies - emin))


def fit_parameters(func, names, x, y, cfg, angles=()):
    """curve_fit over the parameters ``names``.

    cfg["bounds"][name] is [min, max], or a number that fixes the parameter. cfg["p0"][name] is
    the start value (default: the middle of the bounds). Parameters in ``angles`` are given in
    degrees. Returns the values of all parameters, angles in radians.
    """
    bounds, p0 = cfg["bounds"], cfg.get("p0") or {}

    def rad(name, v):
        return np.deg2rad(v) if name in angles else v

    fixed = {n: rad(n, bounds[n]) for n in names if np.isscalar(bounds[n])}
    free = [n for n in names if n not in fixed]
    lo = [rad(n, bounds[n][0]) for n in free]
    hi = [rad(n, bounds[n][1]) for n in free]
    start = [rad(n, p0[n]) if n in p0 else (l + h) / 2 for n, l, h in zip(free, lo, hi)]

    def f(x, *p):
        return func(x, **fixed, **dict(zip(free, p)))

    popt, _ = curve_fit(f, x, y, p0=start, bounds=(lo, hi))
    values = dict(fixed, **dict(zip(free, popt)))
    return [values[n] for n in names]


def report_deviations(de):
    print("Maximal deviation in energy: {:12.6f} (meV)".format(np.max(np.abs(de))))
    print("Root mean square of energy differences: {:12.6f} (meV)\n".format(np.sqrt(np.mean(np.square(de)))))


def cart2sph_phi360(xyz):
    """cart2sph in degrees with phi in [0, 360) and phi = 0 at the poles."""
    r, theta, phi = cart2sph(xyz)
    if abs(np.radians(theta)) < 1e-6 or abs(np.radians(theta) - np.pi) < 1e-6:
        phi = 0.0
    else:
        phi %= 360
        if phi >= 359.95:    # print 0.0 rather than 360.0
            phi = 0.0
    return [r, theta, phi]


# ---------------------------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------------------------

def average(cfg):
    """Average E(n) and E(-n), which removes the exchange energy of a spin rotated against fixed
    neighbours. Writes <data>_average.dat (theta, phi, energy) for directions whose opposite is
    also present."""
    data, (ci, cj, ck) = load(cfg)
    data[data[:, cj] == 360.0, cj] = 0.0
    key = lambda theta, phi: ("{:6.1f}".format(theta), "{:6.1f}".format(phi))
    energies = {key(row[ci], row[cj]): row[ck] for row in data}
    with open(stem(cfg) + "_average.dat", "w") as f:
        for row in data:
            key1 = key(row[ci], row[cj])
            key2 = key(*opposite_direction((row[ci], row[cj])))
            if key1 in energies and key2 in energies:
                f.write("{:9s} {:9s} {:20.12f}\n".format(*key1, (energies[key1] + energies[key2]) / 2))


def change_frame(cfg):
    """Write new_frame.dat: the directions in a frame whose z axis is ``z_direction``, or the
    direction of the lowest (z_direction_from_data: 1) or highest (2) energy."""
    data, (ci, cj, ck) = load(cfg)
    data = data[np.argsort(data[:, ck], kind="stable")]
    z_direction = tuple(cfg.get("z_direction", (0., 0.)))
    mode = cfg.get("z_direction_from_data", 0)
    if mode == 1:
        z_direction = (data[0, ci], data[0, cj])
    elif mode == 2:
        z_direction = (data[-1, ci], data[-1, cj])
    emat_new = get_emat_local(z_direction)
    with open("new_frame.dat", "w") as f:
        for row in data:
            r, theta, phi = change_frame_sph([1, row[ci], row[cj]], np.eye(3), emat_new)
            f.write("{:6.1f} {:6.1f} {:15.9f}\n".format(theta, phi, row[ck]))


def mae(cfg):
    """Out-of-plane and in-plane anisotropy energies from theta = 0 and the theta = 90 circle."""
    data, (ci, cj, ck) = load(cfg)
    pole = data[np.isclose(data[:, ci], 0.0), ck]
    equator = data[np.isclose(data[:, ci], 90.0), ck]
    if len(pole) == 0 or len(equator) == 0:
        raise ValueError(f"mae needs a theta = 0 row and theta = 90 rows in {cfg['data']}.")
    Ez = pole[-1]
    if Ez < np.nanmin(equator):   # easy axis
        mae_out = units.MEV_PER_EV * (np.nanmax(equator) - Ez)
    else:
        mae_out = units.MEV_PER_EV * (Ez - np.nanmin(equator))
    mae_in = units.MEV_PER_EV * (np.nanmax(equator) - np.nanmin(equator))
    print("Out-of-plane MAE: {:6.3f} meV.".format(mae_out))
    print("    In-plane MAE: {:6.3f} meV.".format(mae_in))
    return mae_out, mae_in


def polar_plot(cfg):
    """<data>.pdf: energy (cm^-1, or meV if wavenumber is false) on a (phi, theta) polar map."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    data, (ci, cj, ck) = load(cfg)
    e = relative_meV(data[:, ck], cfg.get("factor", 1.0))
    wavenumber = cfg.get("wavenumber", True)
    if wavenumber:
        e = units.CM_PER_MEV * e
    fig = plt.figure(figsize=(4.2, 3.2))
    ax = fig.add_subplot(projection="polar")
    ax.set_title(cfg.get("title", ""), {"fontsize": 13})
    ax.set_xlim([0, 2 * np.pi])
    ax.set_ylim([0, np.pi])
    ax.set_xticks(np.linspace(0, 2 * np.pi, 9), labels=[r"$\phi=0$", "", r"$\phi=\pi/2$", "", r"$\phi=\pi$", "",
                                                         r"$\phi=3\pi/2$", "", ""])
    ax.set_yticks(np.linspace(0, np.pi, 7), labels=[r"$\theta=0$", "", "", r"$\theta=\pi/2$", "", "", r"$\theta=\pi$"])
    scat = ax.scatter(np.deg2rad(data[:, cj]), np.deg2rad(data[:, ci]), s=60, c=e,
                      vmin=np.nanmin(e), vmax=np.nanmax(e), cmap=mpl.colormaps["rainbow"])
    unit = r"\mathrm{cm}^{-1}" if wavenumber else r"\mathrm{meV}"
    plt.colorbar(scat, label=rf"$E(\theta, \phi)\; ({unit})$", pad=0.15)
    fig.tight_layout()
    fig.savefig(stem(cfg) + ".pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------------------------

def load_for_fit(cfg, select=0):
    """Rows sorted by energy, angles in radians and energies in meV relative to the minimum.
    select = 1 keeps only theta = 0, 90 and 180."""
    data, (ci, cj, ck) = load(cfg)
    if select == 1:
        data = data[np.isin(data[:, ci], (0., 90., 180.))]
    data = data[np.argsort(data[:, ck], kind="stable")]
    e = relative_meV(data[:, ck], cfg.get("factor", 1.0))
    return np.deg2rad(data[:, ci]), np.deg2rad(data[:, cj]), e


def write_fitted(thetas, phis, e, e_fit, fname="energies_fitted.dat"):
    with open(fname, "w") as f:
        for row in zip(np.rad2deg(thetas), np.rad2deg(phis), e, e_fit, e_fit - e):
            f.write("{:6.1f}  {:6.1f}  {:8.3f}   {:8.3f}  {:8.3f}\n".format(*row))


def fit_DE(cfg):
    """E = D [Sz^2 - S(S + 1)/3] + E (Sx^2 - Sy^2) + E0 in the eigenframe; spins as quantum operators (model:
    quantum) or classical vectors (model: classical). Writes energies_fitted.dat,
    d_spin-orbit.dat (diagonal D tensor and the eigenframe axes) and fitting_quality.pdf."""
    frame = np.array(cfg.get("frame", np.eye(3)), dtype=float)
    S = float(cfg["spin"])
    ops = spin_operators(S)
    Sx2, Sy2, Sz2 = ops.Sx @ ops.Sx, ops.Sy @ ops.Sy, ops.Sz @ ops.Sz
    op_D = 3 * Sz2 - ops.ss * ops.ID
    op_E = Sx2 - Sy2

    def quantum(angles, D, E, E0, alpha, beta, gamma):
        thetas, phis = change_frame_angles(angles, frame, rotate_frame(frame, alpha, beta, gamma))
        energies = []
        for theta, phi in zip(thetas, phis):
            n = sph2cart([1, theta, phi], rad=True)
            state = top_state(n[0] * ops.Sx + n[1] * ops.Sy + n[2] * ops.Sz)
            energies.append(D / 3 * expectation(state, op_D) + E * expectation(state, op_E) + E0)
        return np.array(energies)

    def classical(angles, D, E, E0, alpha, beta, gamma):
        thetas, phis = change_frame_angles(angles, frame, rotate_frame(frame, alpha, beta, gamma))
        return S**2 * (D * np.cos(thetas)**2 + E * np.sin(thetas)**2 * np.cos(2 * phis)) + E0

    model = {"quantum": quantum, "classical": classical}[cfg.get("model", "quantum")]
    names = ["D", "E", "E0", "alpha", "beta", "gamma"]
    thetas, phis, e = load_for_fit(cfg, cfg.get("select", 0))
    D, E, E0, alpha, beta, gamma = fit_parameters(model, names, (thetas, phis), e, cfg,
                                                  angles=("alpha", "beta", "gamma"))
    euler = np.rad2deg([alpha, beta, gamma])
    print("D = {:12.6f} (meV), E = {:12.6f} (meV), E0 = {:12.6f} (meV)".format(D, E, E0))
    print("D = {:12.3f} (cm^-1), E = {:12.3f} (cm^-1), E0 = {:12.6f} (cm^-1)".format(
        *(units.CM_PER_MEV * np.array([D, E, E0]))))
    print("Euler angles of the eigenframe: alpha = {:6.1f}, beta = {:6.1f}, gamma = {:6.1f} (deg)\n".format(*euler))

    e_fit = model((thetas, phis), D, E, E0, alpha, beta, gamma)
    report_deviations(e_fit - e)
    write_fitted(thetas, phis, e, e_fit)

    emat = rotate_frame(frame, alpha, beta, gamma)
    print("Eigenframe in spherical coordinates:")
    for axis, v in zip("xyz", emat):
        print("Local e_{}: (theta, phi) = ({:6.1f}, {:6.1f}) deg.".format(axis, *cart2sph_phi360(v)[1:]))
    print("\nEigenframe in cartesian coordinates:")
    print((3 * "{:12.6f} {:12.6f} {:12.6f}\n").format(*emat.flatten()))
    # D tensor = sum_k d_k e_k e_k^T
    d = [-D / 3 + E, -D / 3 - E, 2 * D / 3]
    with open("d_spin-orbit.dat", "w") as f:
        for dk, ek in zip(d, emat):
            f.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(dk, *ek))
    plot_fit(e, e_fit, "fitting_quality.pdf", "magnetic configuration", "E (meV)", cfg.get("title", ""))


def fit_D_tensor(cfg):
    """E = S^T D S + E0 with a general symmetric D (quantum spins), then D and E in the
    eigenframe. Writes energies_fitted.dat and fitting_quality.pdf."""
    ops = spin_operators(float(cfg["spin"]))
    Sx, Sy, Sz = ops.Sx, ops.Sy, ops.Sz
    pair_ops = [Sx @ Sx, Sy @ Sy, Sz @ Sz, Sx @ Sy + Sy @ Sx, Sx @ Sz + Sz @ Sx, Sy @ Sz + Sz @ Sy]

    def model(angles, Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, E0):
        op = sum(Dk * Ok for Dk, Ok in zip((Dxx, Dyy, Dzz, Dxy, Dxz, Dyz), pair_ops))
        energies = []
        for theta, phi in zip(*angles):
            n = sph2cart([1, theta, phi], rad=True)
            energies.append(expectation(top_state(n[0] * Sx + n[1] * Sy + n[2] * Sz), op) + E0)
        return np.array(energies)

    names = ["Dxx", "Dyy", "Dzz", "Dxy", "Dxz", "Dyz", "E0"]
    thetas, phis, e = load_for_fit(cfg)
    p = fit_parameters(model, names, (thetas, phis), e, cfg)
    print("\nThe fitted D tensor matrix elements are:")
    for name, value in zip(names, p):
        print("{:>3s} = {:12.6f} (meV)".format(name, value))
    e_fit = model((thetas, phis), *p)
    print()
    report_deviations(e_fit - e)
    write_fitted(thetas, phis, e, e_fit)

    Dxx, Dyy, Dzz, Dxy, Dxz, Dyz = p[:6]
    w, v = np.linalg.eigh(np.array([[Dxx, Dxy, Dxz], [Dxy, Dyy, Dyz], [Dxz, Dyz, Dzz]]))
    order = np.argsort(w)
    # The axis furthest from the middle eigenvalue is z
    i, j, k = (1, 2, 0) if 4 * w[order[1]] - 2 * w[order[2]] - 2 * w[order[0]] >= 0 else (1, 0, 2)
    D = w[order[k]] - (w[order[i]] + w[order[j]]) / 2
    E = (w[order[i]] - w[order[j]]) / 2
    print("In the eigenframe:")
    print("D = {:9.3f} (meV), E = {:9.3f} (meV)".format(D, E))
    print("D = {:9.3f} (cm^-1), E = {:9.3f} (cm^-1)".format(units.CM_PER_MEV * D, units.CM_PER_MEV * E))
    ex, ey, ez = v[:, order[i]], v[:, order[j]], v[:, order[k]]
    if np.dot(np.cross(ex, ey), ez) < 0:
        ez = -ez
    print("\nAxes of the eigenframe:")
    for axis, vec in zip("xyz", (ex, ey, ez)):
        print("e{} = ( {:12.6f}, {:12.6f}, {:12.6f} )".format(axis, *vec))
    print()
    plot_fit(e, e_fit, "fitting_quality.pdf", "magnetic configuration", "E (meV)", cfg.get("title", ""))


def fit_curve(cfg):
    """E = a (theta - b)^2 + c along one cut (columns: angle, energy). Writes fitting_quality.pdf."""
    data, (ci, ck) = load(cfg)
    x, e = data[:, ci], relative_meV(data[:, ck])

    def model(x, a, b, c):
        return a * (x - b)**2 + c

    a, b, c = fit_parameters(model, ["a", "b", "c"], x, e, cfg)
    print("a = {:12.6f} meV/deg^2, b = {:6.3f} (deg), c = {:6.3f} meV\n".format(a, b, c))
    report_deviations(model(x, a, b, c) - e)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    plot_fit(e, model(xs, a, b, c), "fitting_quality.pdf", "θ (deg)", "E (meV)", cfg.get("title", ""),
             x=x, x_fit=xs)
    emin = units.MEV_PER_EV * np.nanmin(data[:, ck])
    if a > 0:
        print("Fit minimal energy: {:16.3f} meV\n".format(np.min(model(xs, a, b, c)) + emin))
    else:
        print("Fit maximal energy: {:16.3f} meV\n".format(np.max(model(xs, a, b, c)) + emin))


def fit_surface(cfg):
    """E = a theta'^2 + c, theta' measured from the axis (beta, gamma). Writes fitting_quality.pdf."""
    frame = np.array(cfg.get("frame", np.eye(3)), dtype=float)

    def model(angles, a, c, beta, gamma):
        thetas, _ = change_frame_angles(angles, frame, rotate_frame(frame, 0.0, beta, gamma))
        return a * thetas**2 + c

    thetas, phis, e = load_for_fit(cfg)
    a, c, beta, gamma = fit_parameters(model, ["a", "c", "beta", "gamma"], (thetas, phis), e, cfg,
                                       angles=("beta", "gamma"))
    print("a = {:12.6f} meV/rad^2, c = {:12.6f} meV".format(a, c))
    print("Axis: (theta, phi) = (beta, gamma) = ({:6.3f}, {:6.3f}) (deg)\n".format(*np.rad2deg([beta, gamma])))
    e_fit = model((thetas, phis), a, c, beta, gamma)
    report_deviations(e_fit - e)
    plot_fit(e, e_fit, "fitting_quality.pdf", "magnetic configuration", "E (meV)", cfg.get("title", ""))

