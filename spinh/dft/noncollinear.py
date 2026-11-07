"""Noncollinear VASP jobs for sets of spin directions.

A configuration lists a direction (theta, phi) in degrees for each site. Configurations are stored
in the global frame; each site can also have a local frame (rows = axes), and the local
directions are always computed from the current frames. Each configuration runs in a directory
named after its global directions snapped to a grid (see spinh.sphere), e.g. "4319101_0_0".

SiteJobs     one MAGMOM vector per site (and M_CONSTR / RWIGS); INCAR fields {mstring}, {rwigs}.
             directions_of = "magmom" or "spin" (spin = opposite of the moment).
SaxisJobs    one direction for the whole system, written as SAXIS; INCAR field {sstring}. The
             directory names are single integers; WAVECAR is linked from the directory wave_dir.
"""

import os
import re
import subprocess

import numpy as np

from spinh.dft import recipes
from spinh.dft.base import JobSet, action
from spinh.dft.vasp import VaspJob
from spinh.geometry import (angle_between, angle_between_sph, change_frame_sph, get_emat_local,
                            opposite_direction, sph2cart)
from spinh.sphere import Sphere
from spinh.tables import write_direction_table

FRAMES = ("global", "local")


def as_frame(frame):
    """A frame: a 3x3 matrix whose rows are the axes, or a z direction (theta, phi)."""
    frame = np.asarray(frame, dtype=float)
    if frame.shape == (2,):
        return get_emat_local(tuple(frame))
    if frame.shape in ((3, 3), (9,)):
        return frame.reshape(3, 3)
    raise ValueError(f"A frame is a direction (theta, phi) or a 3x3 matrix, not {frame.tolist()}.")


def check_frame_name(frame):
    if frame not in FRAMES:
        raise ValueError(f"frame must be 'global' or 'local', not {frame!r}.")


class DirectionJobs(JobSet):
    """Configurations of site directions, the site frames and the job directories."""

    criterion = "ediff"     # convergence test, see VaspJob.converged
    de0 = 1.e-8

    def __init__(self, n_site, root=None, n_theta=1801, n_phi=3600):
        super().__init__(root)
        self.n_site = n_site
        self.sphere = Sphere(n_theta, n_phi)
        self.frames = [np.eye(3)] * n_site
        self.configurations = []        # global frame
        self.dir_names = []
        self.restart_options = {}

    # -- frames -----------------------------------------------------------------------------------

    def set_frames(self, frames):
        """Local frames: {site: frame} or a list with one frame per site; a frame is a z direction
        (theta, phi) or a 3x3 matrix whose rows are the axes. Unlisted sites use the global frame.
        Stored configurations stay as they are in the global frame."""
        if not isinstance(frames, dict):
            frames = dict(enumerate(frames))
        self.frames = [np.eye(3)] * self.n_site
        for site, frame in frames.items():
            self.frames[site] = as_frame(frame)

    def to_local(self, conf):
        return [change_frame_sph([1, theta, phi], np.eye(3), emat)[1:] for (theta, phi), emat in zip(conf, self.frames)]

    def to_global(self, conf):
        return [change_frame_sph([1, theta, phi], emat, np.eye(3))[1:] for (theta, phi), emat in zip(conf, self.frames)]

    @property
    def configurations_local(self):
        return [self.to_local(conf) for conf in self.configurations]

    # -- configurations ---------------------------------------------------------------------------

    def dir_name(self, conf):
        return "_".join(str(self.sphere.index(d)) for d in conf)

    def add(self, configurations, frame="global"):
        """Add configurations, each a list of (theta, phi), one per site, given in the global or
        the local frames. A configuration whose directory exists already is skipped."""
        check_frame_name(frame)
        for conf in configurations:
            conf = [list(d) for d in conf]
            if len(conf) != self.n_site:
                raise ValueError(f"Configuration {conf} has {len(conf)} directions for {self.n_site} sites.")
            if frame == "local":
                conf = self.to_global(conf)
            name = self.dir_name(conf)
            if name not in self.dir_names:
                self.configurations.append(conf)
                self.dir_names.append(name)

    def add_uniform(self, directions, frame="global"):
        """For each direction, a configuration with every site along it (in its own frame if
        ``frame`` is "local")."""
        self.add([[d] * self.n_site for d in directions], frame)

    def add_dir_name(self, dir_name):
        """Add the configuration of a job directory name."""
        self.add([[list(self.sphere.direction(int(i))) for i in dir_name.split("_")]])

    def set_direction(self, site, direction, frame="global"):
        """Set the direction of one site in every configuration; duplicates are dropped."""
        check_frame_name(frame)
        confs = self.configurations_local if frame == "local" else self.configurations
        for conf in confs:
            conf[site] = list(direction)
        self.configurations, self.dir_names = [], []
        self.add(confs, frame)

    def flip(self, site):
        """Reverse the direction of one site in every configuration."""
        confs = self.configurations
        for conf in confs:
            conf[site] = list(opposite_direction(conf[site]))
        self.configurations, self.dir_names = [], []
        self.add(confs)

    # -- jobs -------------------------------------------------------------------------------------

    def job(self, i):
        return VaspJob(self.root, self.dir_names[i])

    def input_fields(self, conf):
        """INCAR format fields of a configuration (global frame)."""
        raise NotImplementedError

    def wave_dir(self):
        return None

    def write_inputs(self, i):
        self.job(i).write_inputs(wave_dir=self.wave_dir(), **self.input_fields(self.configurations[i]))

    def converged(self, i):
        return self.job(i).converged(self.criterion, self.de0)

    def energy(self, i):
        """Energy of job i (eV), 0 if there is none."""
        return self.job(i).energy()

    def finished(self):
        """Which jobs need no restart."""
        return [self.converged(i)[0] for i in range(self.n_conf)]

    def write_table(self, fname, rows, row_format):
        """fname under root: the global block commented out, then the local block."""
        write_direction_table(self.path(fname), self.configurations,
                              self.configurations_local, rows, row_format)

    # -- actions ----------------------------------------------------------------------------------

    @action("print", "print the configurations (local frames) and their directories")
    def print_configurations(self, frame="local"):
        confs = self.configurations_local if frame == "local" else self.configurations
        fmt = "{:6.1f} {:6.1f} " * self.n_site + ": {:30s}"
        for conf, dir_name in zip(confs, self.dir_names):
            print(fmt.format(*[a for d in conf for a in d], dir_name))

    @action("print-global", "print the configurations in the global frame")
    def print_global(self):
        self.print_configurations("global")

    @action("restart-test", "print what restart would do")
    def restart_test(self):
        self.restart(test=True)

    @action("restart", "restart the unfinished jobs (restart_options)")
    def restart(self, test=False):
        """Resubmit the unfinished jobs. restart_options: with from_neighbor (default True), first
        copy (copy_it, default True) or link file_name (default CHGCAR) from the finished job with
        the closest directions, if all sites are within max_angle (default 180)."""
        options = dict(max_angle=180, from_neighbor=True, file_name="CHGCAR", copy_it=True)
        options.update(self.restart_options)
        finished = self.finished()
        done = [j for j in range(self.n_conf) if finished[j]]
        local = self.configurations_local
        fmt = "({:5.1f},{:5.1f})," * self.n_site
        for i in range(self.n_conf):
            if finished[i]:
                continue
            plan = "Restart (" + fmt.format(*[a for d in local[i] for a in d]) + ") / " + self.dir_names[i]
            if options["from_neighbor"]:
                angles = [max(angle_between_sph([1, *a], [1, *b]) for a, b in
                              zip(self.configurations[i], self.configurations[j])) for j in done]
                if not done or min(angles) > options["max_angle"]:
                    continue
                j = done[int(np.argmin(angles))]
                plan += " from (" + fmt.format(*[a for d in local[j] for a in d]) + ") / " + self.dir_names[j]
            print(plan)
            if not test:
                if options["from_neighbor"]:
                    orig, dest = self.job(j).file(options["file_name"]), self.job(i).file(options["file_name"])
                    subprocess.run(["cp", orig, dest] if options["copy_it"] else ["ln", "-sf", orig, dest])
                self.job(i).submit()

    @action("energies", "write energies.dat (eV)")
    def write_energies(self):
        energies = []
        for i in range(self.n_conf):
            print(self.dir_names[i])
            energies.append(self.energy(i))
        self.write_table("energies.dat", [[e] for e in energies], "{:15.9f}")
        return energies

    @action("all-energies", "add every job directory found under root, then write energies.dat")
    def write_all_energies(self):
        pattern = re.compile(r"\d+" + r"_\d+" * (self.n_site - 1))
        for name in sorted(os.listdir(self.root)):
            if pattern.fullmatch(name) and os.path.isdir(self.path(name)):
                self.add_dir_name(name)
        self.print_configurations()
        return self.write_energies()


class SiteJobs(DirectionJobs):
    """One MAGMOM direction per site."""

    def __init__(self, magnetic_moments, directions_of="magmom", root=None, n_theta=1801, n_phi=3600):
        super().__init__(len(magnetic_moments), root, n_theta, n_phi)
        if directions_of not in ("magmom", "spin"):
            raise ValueError("directions_of must be 'magmom' or 'spin'")
        self.magnetic_moments = list(magnetic_moments)
        self.directions_of = directions_of
        self.exchange_pairs = []        # (i, j, [(dir_name1, dir_name2), ...])

    def add_exchange_pair(self, i, j, rest=recipes.DOWN, axes_i=(recipes.X, recipes.Y, recipes.Z),
                          axes_j=(recipes.X, recipes.Y, recipes.Z), frame="local"):
        """Configurations for the exchange matrix of sites i and j (0-based): site j along b,
        site i along -a and along +a, for a in axes_i and b in axes_j, the other sites along
        ``rest``. Each (-a, +a) pair is one equation of `spinh fit anisotropic-exchange`; the
        "energies" action writes them to exchange_energies.yaml."""
        conf1s, conf2s = [], []
        for b in axes_j:
            for a in axes_i:
                conf = [list(rest)] * self.n_site
                conf[j] = list(b)
                conf1s.append(conf[:i] + [list(opposite_direction(a))] + conf[i + 1:])
                conf2s.append(conf[:i] + [list(a)] + conf[i + 1:])
        self.add(conf1s, frame)
        self.add(conf2s, frame)
        to_global = self.to_global if frame == "local" else (lambda conf: conf)
        names = [(self.dir_name(to_global(c1)), self.dir_name(to_global(c2))) for c1, c2 in zip(conf1s, conf2s)]
        self.exchange_pairs.append((i, j, names))

    def magmom_directions(self, conf):
        if self.directions_of == "spin":
            return [opposite_direction(d) for d in conf]
        return conf

    def mstring(self, conf):
        """MAGMOM vectors, one site per line."""
        return "".join((3 * "{:15.10f}").format(*sph2cart([m, *d])) + "   \\\n"
                       for m, d in zip(self.magnetic_moments, self.magmom_directions(conf)))

    def rwigs(self):
        """Wigner-Seitz radii of the species, from POTCAR."""
        with open(self.path("POTCAR")) as f:
            return "".join(" " + line.split()[5] for line in f if "RWIGS" in line)

    def input_fields(self, conf):
        return {"mstring": self.mstring(conf), "rwigs": self.rwigs()}

    def write_energies(self):
        energies = super().write_energies()
        if self.exchange_pairs:
            self.write_exchange_energies(energies)
        return energies

    def write_exchange_energies(self, energies, fname="exchange_energies.yaml"):
        """The energy pairs of add_exchange_pair, as the `energies:` entries of
        `spinh fit anisotropic-exchange`; pairs with a missing energy are left out."""
        index = {name: k for k, name in enumerate(self.dir_names)}

        def conf(k):
            angles = ", ".join(f"{a:.6f}" for d in self.configurations[k] for a in d)
            return f"[{angles}, {energies[k]:.9f}]"

        lines = ["# Energy pairs for `spinh fit anisotropic-exchange`, written by `python main.py energies`.",
                 "# conf1, conf2: theta, phi of every site (degrees, global frame), then the energy (eV).", ""]
        missing = 0
        for i, j, pairs in self.exchange_pairs:
            for name1, name2 in pairs:
                k1, k2 = index.get(name1), index.get(name2)
                if k1 is None or k2 is None or energies[k1] == 0 or energies[k2] == 0:
                    missing += 1
                    continue
                lines += [f"- pair: [{i + 1}, {j + 1}]", "  frame: global",
                          f"  conf1: {conf(k1)}", f"  conf2: {conf(k2)}", ""]
        with open(self.path(fname), "w") as f:
            f.write("\n".join(lines))
        print(f"Wrote {self.path(fname)}" + (f"; {missing} pair(s) without energies left out." if missing else "."))

    @action("dipoles", "write dipoles.dat (dipole correction)")
    def write_dipoles(self, crystal=False):
        """dipoles.dat: px, py, pz, |p| (e Ang)."""
        dipoles = []
        for i in range(self.n_conf):
            print(self.dir_names[i])
            dipoles.append(self.job(i).dipole(crystal))
        self.write_table("dipoles.dat", dipoles, "{:12.6f} {:12.6f} {:12.6f} {:12.6f}")
        return dipoles

    @action("dipoles-crystal", "write dipoles.dat with the Berry-phase dipoles")
    def write_dipoles_crystal(self):
        return self.write_dipoles(crystal=True)

    @action("magmoms", "print the angles between the computed and the initial site moments")
    def print_magmom_angles(self):
        for i, conf in enumerate(self.configurations):
            moments = self.job(i).local_moments(self.n_site)
            angles = [np.nan if np.isnan(m).any() else angle_between(m, sph2cart([1, *d]))
                      for m, d in zip(moments, self.magmom_directions(conf))]
            print(("{:50s} " + "{:6.1f} " * self.n_site).format(self.dir_names[i], *angles))

    @action("site-moments", "write lmms_and_loms_site<i>.dat: spin and orbital moments of site i (arg: i)")
    def write_site_moments(self, i_site):
        """Spin and orbital moments of one site and the angle between them."""
        i_site = int(i_site)
        rows = []
        for i in range(self.n_conf):
            print(self.dir_names[i])
            lmm = self.job(i).local_moments(self.n_site, "magnetization")[i_site]
            lom = self.job(i).local_moments(self.n_site, "orbital moment")[i_site]
            angle = np.nan if np.isnan(lmm).any() or np.isnan(lom).any() else angle_between(lmm, lom)
            rows.append([*lmm, *lom, angle])
        self.write_table(f"lmms_and_loms_site{i_site}.dat", rows, " ".join(["{:12.6f}"] * 7))

    @action("occmat", "write occmat and occmat.eigen in every job directory")
    def write_occupation_matrices(self):
        for i in range(self.n_conf):
            self.job(i).occupation_matrices()


class SaxisJobs(DirectionJobs):
    """One spin-quantization axis (SAXIS) for the whole system.

    A job counts as converged when the energy change in its ``output`` is below de0, and its
    energy only when it is also within max_energy of e_ref (otherwise 0 in energies.dat).
    """

    criterion = "output"

    def __init__(self, wave_dir="0", e_ref=0.0, max_energy=np.inf, de0=1.e-8, root=None, n_theta=1801, n_phi=3600):
        super().__init__(1, root, n_theta, n_phi)
        self.wave_dir_name = wave_dir
        self.e_ref = e_ref
        self.max_energy = max_energy
        self.de0 = de0
        self.restart_options = {"from_neighbor": False}

    def set_frame(self, frame):
        """The local frame: a z direction (theta, phi) or a 3x3 matrix whose rows are the axes."""
        self.set_frames([frame])

    def wave_dir(self):
        return self.wave_dir_name

    def input_fields(self, conf):
        return {"sstring": "{:10.6f} {:10.6f} {:10.6f}".format(*sph2cart([1, *conf[0]]))}

    def energy(self, i):
        job = self.job(i)
        e = job.energy()
        if job.converged(self.criterion, self.de0)[0] and abs(e - self.e_ref) < self.max_energy:
            return e
        return 0

    def finished(self):
        return [self.energy(i) != 0 for i in range(self.n_conf)]
