"""Noncollinear VASP jobs for sets of spin directions.

A configuration lists a direction (theta, phi) in degrees for each site. Directions can be given
in the global frame or in local frames (one per site, rows = axes), and are kept in both. Each
configuration runs in a directory named after its directions snapped to a grid (see
spinh.sphere), e.g. "4319101_0_0".

SiteJobs     one MAGMOM vector per site (and M_CONSTR / RWIGS); INCAR fields {mstring}, {rwigs}.
             directions_of = "magmom" or "spin" (spin = opposite of the moment).
SaxisJobs    one direction for the whole system, written as SAXIS; INCAR field {sstring}. The
             directory names are single integers; WAVECAR is linked from data.wave_dir.
"""

import os
import re
import subprocess

import numpy as np

from spinh.dft.slurm import run_action
from spinh.dft.vasp import VaspJob
from spinh.geometry import (angle_between, angle_between_sph, cart2sph, change_frame_sph, get_emat_local,
                          opposite_direction, sph2cart)
from spinh.sphere import Sphere
from spinh.tables import write_direction_table


class DirectionJobs:
    """Configurations of site directions, their frames and job directories."""

    criterion = "ediff"     # convergence test, see VaspJob.converged

    def __init__(self, data, n_site, n_theta=1801, n_phi=3600):
        self.data = data
        self.n_site = n_site
        self.sphere = Sphere(n_theta, n_phi)
        self.z_directions = [(0, 0)] * n_site
        self.emats = [np.eye(3)] * n_site
        self.configurations = []
        self.configurations_local = []
        self.dir_names = []
        self.energies = []
        self.restart_options = {}

    @property
    def n_conf(self):
        return len(self.configurations)

    # -- frames and directory names -------------------------------------------------------------

    def to_local(self, conf):
        return [change_frame_sph([1, theta, phi], np.eye(3), emat)[1:] for (theta, phi), emat in zip(conf, self.emats)]

    def to_global(self, conf_local):
        return [change_frame_sph([1, theta, phi], emat, np.eye(3))[1:] for (theta, phi), emat in zip(conf_local, self.emats)]

    def both_frames(self, conf, local_ref_frame):
        """(global, local) versions of a configuration given in either frame."""
        conf = [list(d) for d in conf]
        return (self.to_global(conf), conf) if local_ref_frame else (conf, self.to_local(conf))

    def dir_name(self, conf):
        return "_".join(str(self.sphere.index(d)) for d in conf)

    def conf_of_dir_name(self, dir_name):
        return [list(self.sphere.direction(int(i))) for i in dir_name.split("_")[:self.n_site]]

    def set_frames(self, sites, z_directions=(), emats=(), from_file=False, emats_file=()):
        """Frames of the listed sites (the others: global); does not transform the stored
        configurations, call transform_configurations for that."""
        self.z_directions = [(0, 0)] * self.n_site
        self.emats = [np.eye(3)] * self.n_site
        for k, site in enumerate(sites):
            if from_file:
                self.emats[site] = np.array(emats_file[site], dtype=float)
            elif len(emats):
                self.emats[site] = np.array(emats[k], dtype=float)
            else:
                self.z_directions[site] = z_directions[k]
                self.emats[site] = get_emat_local(z_direction=z_directions[k])

    # -- configurations ---------------------------------------------------------------------------

    def add_configurations(self, configurations=(), local_ref_frame=False):
        """Add configurations whose directory is new."""
        for conf in configurations:
            conf_global, conf_local = self.both_frames(conf, local_ref_frame)
            dir_name = self.dir_name(conf_global)
            if dir_name not in self.dir_names:
                self.configurations.append(conf_global)
                self.configurations_local.append(conf_local)
                self.dir_names.append(dir_name)
                self.energies.append(0)

    def add_one_configuration_by_dir_name(self, dir_name):
        self.add_configurations([self.conf_of_dir_name(dir_name)], local_ref_frame=False)

    def add_thetas(self, phi=0, theta_min=0, theta_max=180, ntheta=2, local_ref_frame=False):
        """All sites along (theta, phi) for ``ntheta`` thetas at fixed phi."""
        for theta in np.linspace(theta_min, theta_max, ntheta, endpoint=True):
            self.add_configurations([[[theta, phi]] * self.n_site], local_ref_frame)

    def add_phis(self, theta=0, phi_min=0, phi_max=360, nphi=2, local_ref_frame=False):
        """All sites along (theta, phi) for ``nphi`` phis at fixed theta."""
        for phi in np.linspace(phi_min, phi_max, nphi, endpoint=True):
            self.add_configurations([[[theta, phi]] * self.n_site], local_ref_frame)

    def transform_configurations(self, global2local=True):
        """Recompute the local directions from the global ones (or the reverse) with the current
        frames."""
        for i in range(self.n_conf):
            if global2local:
                self.configurations_local[i] = self.to_local(self.configurations[i])
            else:
                self.configurations[i] = self.to_global(self.configurations_local[i])
        self.dir_names = [self.dir_name(conf) for conf in self.configurations]

    def set_direction_for_one_site(self, i_site=0, direction=(0, 0), local_ref_frame=False):
        """Set the direction of one site in every configuration."""
        self.transform_configurations(global2local=True)
        confs = self.configurations_local if local_ref_frame else self.configurations
        for conf in confs:
            conf[i_site] = list(direction)
        self.transform_configurations(global2local=not local_ref_frame)

    def flip_direction_for_one_site(self, i_site=0):
        """Reverse the direction of one site in every configuration."""
        for conf in self.configurations:
            conf[i_site] = list(opposite_direction(conf[i_site]))
        self.transform_configurations(global2local=True)

    def print_dirs_and_configs(self, local_ref_frame=False):
        confs = self.configurations_local if local_ref_frame else self.configurations
        fmt = "{:6.1f} {:6.1f} " * self.n_site + ": {:30s}"
        for conf, dir_name in zip(confs, self.dir_names):
            print(fmt.format(*[a for d in conf for a in d], dir_name))

    # -- jobs -------------------------------------------------------------------------------------

    def job(self, i):
        return VaspJob(self.data, self.dir_names[i])

    def input_fields(self, conf):
        """INCAR format fields of a configuration (global frame)."""
        raise NotImplementedError

    def wave_dir(self):
        return None

    def setup_jobs(self, submit=False):
        """Create the directories and write the inputs (or submit, if ``submit``)."""
        for i, conf in enumerate(self.configurations):
            job = self.job(i)
            job.create_dir()
            if submit:
                job.submit()
            else:
                job.write_inputs(wave_dir=self.wave_dir(), **self.input_fields(conf))

    def submit_jobs(self):
        self.setup_jobs(submit=True)

    def check_convergences(self, restart=False, de0=1.e-8):
        self.convergences = [self.job(i).check_convergence(restart, self.criterion, de0) for i in range(self.n_conf)]
        return self.convergences

    def done(self):
        """Which jobs need no restart."""
        return self.check_convergences()

    def restart(self, test=True, max_angle=180, from_neighbor=True, file_name="CHGCAR", copy_it=True):
        """Resubmit the unfinished jobs. With ``from_neighbor``, first copy (or link) ``file_name``
        from the finished job with the closest directions, if all sites are within ``max_angle``.
        ``test`` only prints the plan."""
        done = self.done()
        finished = [j for j in range(self.n_conf) if done[j]]
        fmt = "({:5.1f},{:5.1f})," * self.n_site
        for i in range(self.n_conf):
            if done[i]:
                continue
            plan = "Restart (" + fmt.format(*[a for d in self.configurations_local[i] for a in d]) + ") / " + self.dir_names[i]
            if from_neighbor:
                angles = [max(angle_between_sph([1, *a], [1, *b]) for a, b in
                              zip(self.configurations[i], self.configurations[j])) for j in finished]
                if not finished or min(angles) > max_angle:
                    continue
                j = finished[int(np.argmin(angles))]
                plan += " from (" + fmt.format(*[a for d in self.configurations_local[j] for a in d]) + ") / " + self.dir_names[j]
            print(plan)
            if not test:
                if from_neighbor:
                    orig, dest = self.job(j).file(file_name), self.job(i).file(file_name)
                    subprocess.run(["cp", orig, dest] if copy_it else ["ln", "-sf", orig, dest])
                self.job(i).submit()

    def energy(self, i):
        return self.job(i).energy()

    def get_energies(self):
        """Write energies.dat (eV): the global frame commented out, then the local frame."""
        self.energies = []
        for i in range(self.n_conf):
            print(self.dir_names[i])
            self.energies.append(self.energy(i))
        write_direction_table(os.path.join(self.data.root_dir, "energies.dat"), self.configurations,
                              self.configurations_local, [[e] for e in self.energies], "{:15.9f}")

    def get_all_energies(self):
        """Add every job directory under root_dir, then get_energies."""
        pattern = re.compile(r"\d+" + r"_\d+" * (self.n_site - 1))
        for name in sorted(os.listdir(self.data.root_dir)):
            if pattern.fullmatch(name) and os.path.isdir(os.path.join(self.data.root_dir, name)):
                self.add_one_configuration_by_dir_name(name)
        self.print_dirs_and_configs(local_ref_frame=True)
        self.get_energies()

    def actions(self):
        """Actions of ``python main.py ACTION``."""
        return {
            "print": lambda: self.print_dirs_and_configs(local_ref_frame=True),
            "print-global": lambda: self.print_dirs_and_configs(local_ref_frame=False),
            "setup": self.setup_jobs,
            "submit": self.submit_jobs,
            "check": self.check_convergences,
            "resubmit": lambda: self.check_convergences(restart=True),
            "restart-test": lambda: self.restart(test=True, **self.restart_options),
            "restart": lambda: self.restart(test=False, **self.restart_options),
            "energies": self.get_energies,
            "all-energies": self.get_all_energies,
        }

    def run(self, argv=None):
        """Run the action named on the command line (default: print)."""
        return run_action(self.actions(), argv, ACTION_DOCS)


ACTION_DOCS = {
    "print": "print the configurations (local frames) and their directories",
    "print-global": "print the configurations in the global frame",
    "setup": "create the directories and write the inputs",
    "submit": "submit the jobs",
    "check": "print the convergence of every job",
    "resubmit": "resubmit the jobs that are not converged",
    "restart-test": "print what restart would do",
    "restart": "restart unfinished jobs (restart_options)",
    "energies": "write energies.dat",
    "all-energies": "write energies.dat for every job directory found",
    "dipoles": "write dipoles.dat",
    "dipoles-crystal": "write dipoles.dat with the Berry-phase dipoles",
    "check-magmoms": "print the angles between the computed and the initial moments",
    "moments": "write lmms_and_loms_site<i>.dat (args: site index)",
    "occmat": "write occmat and occmat.eigen in every job directory",
}


class SiteJobs(DirectionJobs):
    """One MAGMOM direction per site."""

    def __init__(self, data, n_site=3, magnetic_moments=(3, 3, 3), directions_of="magmom",
                 n_theta=1801, n_phi=3600):
        super().__init__(data, n_site, n_theta, n_phi)
        if directions_of not in ("magmom", "spin"):
            raise ValueError("directions_of must be 'magmom' or 'spin'")
        self.magnetic_moments = list(magnetic_moments)
        self.directions_of = directions_of
        self.dipoles = []

    # Names used by existing main.py files
    add_thetas_colinear_magmom = add_thetas_colinear_spin = DirectionJobs.add_thetas
    add_phis_colinear_magmom = add_phis_colinear_spin = DirectionJobs.add_phis

    def setup_local_ref_frames(self, sites=(), z_directions=(), from_file=False, emats=(), from_z_direction=True):
        """Local frames of the listed sites: along z_directions, from data.emats_file, or
        ``emats`` (with from_z_direction=False)."""
        self.set_frames(sites, z_directions, () if from_z_direction else emats, from_file,
                        getattr(self.data, "emats_file", ()))

    def magmom_directions(self, conf):
        if self.directions_of == "spin":
            return [opposite_direction(d) for d in conf]
        return conf

    def mstring(self, conf):
        """MAGMOM vectors, one site per line."""
        return "".join((3 * "{:15.10f}").format(*sph2cart([m, *d])) + "   \\\n"
                       for m, d in zip(self.magnetic_moments, self.magmom_directions(conf)))

    def input_fields(self, conf):
        return {"mstring": self.mstring(conf), "rwigs": self.rwigs()}

    def rwigs(self):
        """Wigner-Seitz radii of the species, from POTCAR."""
        with open(os.path.join(self.data.root_dir, "POTCAR")) as f:
            return "".join(" " + line.split()[5] for line in f if "RWIGS" in line)

    def get_dipoles(self, crystal=False):
        """Write dipoles.dat (px, py, pz, |p|; e Ang)."""
        self.dipoles = []
        for i in range(self.n_conf):
            print(self.dir_names[i])
            self.dipoles.append(self.job(i).dipole(crystal))
        write_direction_table(os.path.join(self.data.root_dir, "dipoles.dat"), self.configurations,
                              self.configurations_local, self.dipoles, "{:12.6f} {:12.6f} {:12.6f} {:12.6f}")

    def check_local_magmoms_all_configurations(self):
        """Print, per job, the angles between the computed and the initial site moments."""
        for i, conf in enumerate(self.configurations):
            moments = self.job(i).local_moments(self.n_site)
            angles = [np.nan if np.isnan(m).any() else angle_between(m, sph2cart([1, *d]))
                      for m, d in zip(moments, self.magmom_directions(conf))]
            print(("{:50s} " + "{:6.1f} " * self.n_site).format(self.dir_names[i], *angles))

    def get_local_magmom_and_orbmom_all_configurations(self, i_site):
        """Write lmms_and_loms_site<i>.dat: spin and orbital moments of one site and their angle."""
        i_site = int(i_site)
        rows = []
        for i in range(self.n_conf):
            print(self.dir_names[i])
            lmm = self.job(i).local_moments(self.n_site, "magnetization")[i_site]
            lom = self.job(i).local_moments(self.n_site, "orbital moment")[i_site]
            angle = np.nan if np.isnan(lmm).any() or np.isnan(lom).any() else angle_between(lmm, lom)
            rows.append([*lmm, *lom, angle])
        write_direction_table(os.path.join(self.data.root_dir, f"lmms_and_loms_site{i_site}.dat"),
                              self.configurations, self.configurations_local, rows, " ".join(["{:12.6f}"] * 7))

    def get_occmat_all_configurations(self):
        for i in range(self.n_conf):
            self.job(i).occupation_matrices()

    def actions(self):
        actions = super().actions()
        actions.update({
            "dipoles": self.get_dipoles,
            "dipoles-crystal": lambda: self.get_dipoles(crystal=True),
            "check-magmoms": self.check_local_magmoms_all_configurations,
            "moments": self.get_local_magmom_and_orbmom_all_configurations,
            "occmat": self.get_occmat_all_configurations,
        })
        return actions


class SaxisJobs(DirectionJobs):
    """One spin-quantization axis (SAXIS) for the whole system."""

    criterion = "output"

    def __init__(self, data, n_theta=1801, n_phi=3600):
        super().__init__(data, 1, n_theta, n_phi)
        self.e_ref = 0.0          # energies further than max_energy from e_ref count as failed
        self.max_energy = np.inf
        self.de0 = 1.e-8
        self.restart_options = {"from_neighbor": False}

    @property
    def directions(self):
        return [conf[0] for conf in self.configurations]

    @property
    def directions_local(self):
        return [conf[0] for conf in self.configurations_local]

    def setup_local_ref_frames(self, z_direction=(0, 0), from_file=False):
        """Local frame along ``z_direction``, or data.emat_file."""
        if from_file:
            self.set_frames([0], emats=[self.data.emat_file])
            self.z_directions = [tuple(cart2sph(self.data.emat_file[2])[1:])]
        else:
            self.set_frames([0], z_directions=[z_direction])

    def add_directions(self, directions=((0, 0),), local_ref_frame=False):
        self.add_configurations([[d] for d in directions], local_ref_frame)

    def add_directions_by_dir_name(self, dir_name):
        self.add_one_configuration_by_dir_name(dir_name)

    def transform_directions(self, global2local=True):
        self.transform_configurations(global2local)

    def wave_dir(self):
        return self.data.wave_dir

    def input_fields(self, conf):
        return {"sstring": "{:10.6f} {:10.6f} {:10.6f}".format(*sph2cart([1, *conf[0]]))}

    def check_convergences(self, restart=False, de0=None):
        return super().check_convergences(restart, self.de0 if de0 is None else de0)

    def energy(self, i):
        """The energy if converged and within max_energy of e_ref, otherwise 0."""
        job = self.job(i)
        e = job.energy()
        if job.converged(self.criterion, self.de0)[0] and abs(e - self.e_ref) < self.max_energy:
            return e
        return 0

    def done(self):
        self.get_energies()
        return [e != 0 for e in self.energies]
