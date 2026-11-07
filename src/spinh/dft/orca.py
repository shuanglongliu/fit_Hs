"""Broken-symmetry ORCA jobs, one directory per 'u'/'d' configuration.

``data`` (the job directory's data.py) provides root_dir, hs_dir (high-spin run whose orca.gbw is
the starting guess), spin_only_Szs, effective_spins, expected_magmoms, configs, n_cpus,
mem_per_core, input_file_hs, input_file_bs and job_script. The input templates have the fields
{n_cpus}, {mem_per_core}, {flip_spin_tag}, {flip_spin}, {final_ms_tag}, {final_ms} and
{multiplicity}; orca.xyz is linked from root_dir.

Workflow: setup, submit, check, then for unconverged runs backup (archives run_N.* and points
dft.gbw at the last orbitals), restart-inputs (input without FlipSpin/FinalMs) and resubmit;
finally magmoms and energies.
"""

import csv
import os
import subprocess

from spinh.dft.slurm import run_action
from spinh.tables import write_collinear_energies


class BsDft:
    """One broken-symmetry calculation."""

    def __init__(self, data, config="ud", restart=False):
        self.data = data
        self.config = config
        self.spin_only_Szs = data.spin_only_Szs
        self.restart = restart
        self.work_dir = os.path.join(data.root_dir, config)
        self.n_sites = len(config)
        self.high_spin = all(c == "u" for c in config)
        self.input_file = data.input_file_hs if self.high_spin else data.input_file_bs
        self.energy = None
        self.correct_magmoms = False
        self._get_settings()
        self.ref_gbw_file = self._reference_gbw()

    def _get_settings(self):
        """FlipSpin / FinalMs directives and the multiplicity."""
        # A restart reads the broken-symmetry state from the GBW file, so no spin flip
        prefix = "# " if self.restart else ""
        self.flip_spin_tag, self.final_ms_tag = prefix + "FlipSpin", prefix + "FinalMs"
        # The computed state has all spins reversed if there are more down than up spins, so that
        # the total Ms is positive; the directory keeps the name of the configuration.
        self.state = self.config
        if self.state.count("d") > self.state.count("u"):
            self.state = self.state.replace("u", "x").replace("d", "u").replace("x", "d")
        self.flip_spin = ", ".join(str(i) for i, c in enumerate(self.state) if c == "d")
        self.final_ms = sum(s if c == "u" else -s for c, s in zip(self.state, self.spin_only_Szs))
        S = abs(self.final_ms) if self.restart else sum(self.spin_only_Szs)
        self.multiplicity = int(2 * S + 1)

    def _run_files(self, ext):
        if not os.path.isdir(self.work_dir):
            return []
        return sorted(int(f[4:-len(ext)]) for f in os.listdir(self.work_dir)
                      if f.startswith("run_") and f.endswith(ext))

    def _reference_gbw(self):
        """Starting orbitals: the high-spin run, or for a restart the latest run_N.gbw."""
        runs = self._run_files(".gbw") if self.restart else []
        if runs:
            return os.path.join(self.work_dir, f"run_{runs[-1]}.gbw")
        return os.path.join(self.data.hs_dir, "orca.gbw")

    def create_work_dir(self):
        os.makedirs(self.work_dir, exist_ok=True)

    def get_input_file(self):
        """Write orca.inp and orca.job; link orca.xyz and (except high spin) dft.gbw."""
        print(self.work_dir)
        with open(os.path.join(self.work_dir, "orca.inp"), "w") as f:
            f.write(self.input_file.format(
                n_cpus=self.data.n_cpus, mem_per_core=self.data.mem_per_core,
                flip_spin_tag=self.flip_spin_tag, flip_spin=self.flip_spin,
                final_ms_tag=self.final_ms_tag, final_ms=self.final_ms, multiplicity=self.multiplicity))
        with open(os.path.join(self.work_dir, "orca.job"), "w") as f:
            f.write(self.data.job_script)
        subprocess.run(["ln", "-sf", "../orca.xyz", os.path.join(self.work_dir, "orca.xyz")])
        if not self.high_spin:
            subprocess.run(["ln", "-sf", self.ref_gbw_file, os.path.join(self.work_dir, "dft.gbw")])

    def submit_job(self):
        print(self.work_dir)
        subprocess.run(["sbatch", "orca.job"], cwd=self.work_dir)

    def output(self):
        path = os.path.join(self.work_dir, "orca.out")
        if not os.path.exists(path):
            return None
        with open(path) as f:
            return f.readlines()

    def check_convergence(self, verbose=True):
        """Converged when "SCF CONVERGED" appears twice after a spin flip (high spin, then
        broken symmetry) and once otherwise."""
        lines = self.output()
        self.converged = False
        if lines is not None:
            expected = 2 if "F L I P    S P I N" in "".join(lines) else 1
            self.converged = sum("SCF CONVERGED" in line for line in lines) == expected
        if verbose:
            print(f"Configuration {self.config} is {'' if self.converged else 'not '}converged.")
        return self.converged

    def backup_files(self):
        """For an unconverged run: delete scratch files, archive orca.inp/out/gbw as run_N.*, and
        link dft.gbw to run_N.gbw."""
        if self.check_convergence(verbose=False):
            print(f"Configuration {self.config}: converged, skipping backup ...")
            return
        wd = self.work_dir
        for fname in os.listdir(wd):
            if ".tmp" in fname or ".bas" in fname:
                os.remove(os.path.join(wd, fname))
        gbw_path = os.path.join(wd, "orca.gbw")
        if not os.path.exists(gbw_path):
            return
        n = max(self._run_files(".gbw"), default=0) + 1
        subprocess.run(["cp", os.path.join(wd, "orca.inp"), os.path.join(wd, f"run_{n}.inp")])
        subprocess.run(["cp", os.path.join(wd, "orca.out"), os.path.join(wd, f"run_{n}.out")])
        os.rename(gbw_path, os.path.join(wd, f"run_{n}.gbw"))
        subprocess.run(["ln", "-sf", f"run_{n}.gbw", os.path.join(wd, "dft.gbw")])
        print(f"Configuration {self.config}: backed up the .inp, .out, and .gbw files to run_{n}.* "
              f"and updated dft.gbw symlink.")

    def resubmit_job(self):
        if self.check_convergence(verbose=False):
            print(f"Configuration {self.config}: converged, skipping resubmission ...")
        else:
            print(f"Configuration {self.config}: not converged, resubmitting job ...")
            subprocess.run(["sbatch", "orca.job"], cwd=self.work_dir)

    def check_magnetic_moments(self, moments, magmom_tolerance=0.3):
        """Sign ('u' in the computed state positive) and size (expected_magmoms) of every site's
        spin population."""
        self.correct_magmoms = True
        for moment in moments:
            i = moment["atomic_index"]
            population = moment["spin_population"]
            if population * (1 if self.state[i] == "u" else -1) < 0:
                self.correct_magmoms = False
                print(f"Configuration {self.config}: site {i} (0-based, {moment['element']}) has unexpected "
                      f"spin population in the sign {population}.")
            if abs(abs(population) - self.data.expected_magmoms[i]) > magmom_tolerance:
                self.correct_magmoms = False
                print(f"Configuration {self.config}: site {i} (0-based, {moment['element']}) has unexpected "
                      f"spin population in the size {population} (tolerance: {magmom_tolerance} μB).")

    def get_magnetic_moments(self, magmom_tolerance=0.3):
        """Loewdin spin populations of the magnetic sites from the last section of orca.out;
        checked, then written to magnetic_moments.csv."""
        if not self.check_convergence(verbose=False):
            print(f"Configuration {self.config} is not converged. Stopping ...")
            return
        moments = []
        in_section = False
        for line in self.output():
            if "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" in line:
                moments, in_section = [], True     # keep only the last section
                continue
            if in_section:
                if line.strip() == "" or line.startswith("----"):
                    continue
                if len(moments) < self.n_sites:
                    parts = line.split()
                    moments.append({"atomic_index": int(parts[0]), "element": parts[1].rstrip(":"),
                                    "spin_population": round(float(parts[3]), 2)})
                else:
                    in_section = False
        self.check_magnetic_moments(moments, magmom_tolerance)
        with open(os.path.join(self.work_dir, "magnetic_moments.csv"), "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["atomic_index", "element", "spin_population"])
            writer.writeheader()
            writer.writerows(moments)

    def get_energy(self, magmom_tolerance=0.3):
        """Final single point energy (Ha) if converged with the expected moments, else None."""
        self.energy = None
        if not self.check_convergence(verbose=False):
            return None
        self.get_magnetic_moments(magmom_tolerance)
        if not self.correct_magmoms:
            return None
        for line in reversed(self.output()):
            if "FINAL SINGLE POINT ENERGY" in line:
                self.energy = float(line.split()[-1])
                break
        return self.energy


class BsDftManager:
    """All configurations of data.configs."""

    def __init__(self, data, configs=None, restart=False):
        self.data = data
        self.configs = data.configs if configs is None else configs
        self.jobs = [BsDft(data, config, restart) for config in self.configs]

    def create_work_dirs(self):
        for job in self.jobs:
            job.create_work_dir()

    def get_input_files(self):
        for job in self.jobs:
            job.get_input_file()

    def setup(self):
        self.create_work_dirs()
        self.get_input_files()

    def submit_jobs(self):
        for job in self.jobs:
            job.submit_job()

    def check_convergence(self):
        for job in self.jobs:
            job.check_convergence()

    def backup_files(self):
        for job in self.jobs:
            job.backup_files()

    def resubmit_jobs(self):
        for job in self.jobs:
            job.resubmit_job()

    def get_magnetic_moments(self, magmom_tolerance=0.3):
        for job in self.jobs:
            job.get_magnetic_moments(magmom_tolerance)

    def get_energies(self, magmom_tolerance=0.3):
        """Write energies.dat (Ha) for the isotropic-exchange fit."""
        fname = os.path.join(self.data.root_dir, "energies.dat")
        energies = [job.get_energy(magmom_tolerance) for job in self.jobs]
        write_collinear_energies(fname, [job.state for job in self.jobs], self.data.effective_spins, energies, "Ha")
        print(f"Saved energies for all configurations to {fname}")

    def run(self, argv=None, magmom_tolerance=0.3):
        actions = {
            "setup": self.setup, "submit": self.submit_jobs, "check": self.check_convergence,
            "backup": self.backup_files,
            "restart-inputs": lambda: BsDftManager(self.data, self.configs, restart=True).get_input_files(),
            "resubmit": self.resubmit_jobs,
            "magmoms": lambda tol=magmom_tolerance: self.get_magnetic_moments(float(tol)),
            "energies": lambda tol=magmom_tolerance: self.get_energies(float(tol)),
        }
        docs = {"setup": "create the directories and write the inputs", "submit": "submit the jobs",
                "check": "print the convergence of every job",
                "backup": "archive unconverged runs as run_N.* and link dft.gbw to their orbitals",
                "restart-inputs": "rewrite the inputs for a restart (no FlipSpin/FinalMs)",
                "resubmit": "resubmit the unconverged jobs",
                "magmoms": "check the spin populations [tolerance]",
                "energies": "write energies.dat [tolerance]"}
        return run_action(actions, argv, docs, default="check")
