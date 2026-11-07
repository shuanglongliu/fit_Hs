"""Broken-symmetry ORCA jobs, one directory per 'u'/'d' configuration.

The root directory (that of main.py) holds orca_hs.inp (input of the high-spin run), orca_bs.inp
(input of the broken-symmetry runs), orca.job (the Slurm script) and orca.xyz (linked into every
run). The inputs are templates with the fields {n_cpus}, {mem_per_core}, {flip_spin_tag},
{flip_spin}, {final_ms_tag}, {final_ms} and {multiplicity}. The broken-symmetry runs start from
the orbitals of the high-spin run (hs_dir/orca.gbw).

Workflow: setup, submit, check, then for unconverged runs backup (archives run_N.* and points
dft.gbw at the last orbitals), restart-inputs (input without FlipSpin/FinalMs) and resubmit;
finally magmoms and energies.
"""

import csv
import os
import subprocess

from spinh.dft.base import action
from spinh.dft.collinear import CollinearJobs
from spinh.dft.slurm import sbatch
from spinh.dft.vasp import read_root_file


def spin_state(config):
    """The state ORCA computes: all spins reversed if there are more down than up spins, so that
    the total Ms is positive. The directory keeps the name of the configuration."""
    if config.count("d") > config.count("u"):
        return config.translate(str.maketrans("ud", "du"))
    return config


class OrcaJob:
    """One broken-symmetry calculation in root/<config>; ``jobs`` is its CollinearOrcaJobs."""

    def __init__(self, jobs, config, restart=False):
        self.jobs = jobs
        self.config = config
        self.restart = restart
        self.path = jobs.path(config)
        self.high_spin = all(c == "u" for c in config)
        self.state = spin_state(config)

    def file(self, name):
        return os.path.join(self.path, name)

    def fields(self):
        """FlipSpin / FinalMs directives and the multiplicity. A restart reads the broken-symmetry
        state from the GBW file, so it has no spin flip."""
        prefix = "# " if self.restart else ""
        spins = self.jobs.spins
        final_ms = sum(s if c == "u" else -s for c, s in zip(self.state, spins))
        S = abs(final_ms) if self.restart else sum(spins)
        return dict(n_cpus=self.jobs.n_cpus, mem_per_core=self.jobs.mem_per_core,
                    flip_spin_tag=prefix + "FlipSpin",
                    flip_spin=", ".join(str(i) for i, c in enumerate(self.state) if c == "d"),
                    final_ms_tag=prefix + "FinalMs", final_ms=final_ms, multiplicity=int(2 * S + 1))

    def runs(self):
        """Numbers N of the archived run_N.gbw files."""
        if not os.path.isdir(self.path):
            return []
        return sorted(int(f[4:-4]) for f in os.listdir(self.path) if f.startswith("run_") and f.endswith(".gbw"))

    def reference_gbw(self):
        """Starting orbitals: the high-spin run, or for a restart the latest run_N.gbw."""
        runs = self.runs() if self.restart else []
        if runs:
            return self.file(f"run_{runs[-1]}.gbw")
        return os.path.join(self.jobs.path(self.jobs.hs_dir), "orca.gbw")

    def create_dir(self):
        os.makedirs(self.path, exist_ok=True)

    def write_inputs(self):
        """orca.inp and orca.job; links to orca.xyz and (except high spin) dft.gbw."""
        print(self.path)
        template = read_root_file(self.jobs.root, "orca_hs.inp" if self.high_spin else "orca_bs.inp")
        with open(self.file("orca.inp"), "w") as f:
            f.write(template.format(**self.fields()))
        with open(self.file("orca.job"), "w") as f:
            f.write(read_root_file(self.jobs.root, "orca.job"))
        subprocess.run(["ln", "-sf", "../orca.xyz", self.file("orca.xyz")])
        if not self.high_spin:
            subprocess.run(["ln", "-sf", self.reference_gbw(), self.file("dft.gbw")])

    def submit(self):
        sbatch(self.path, "orca.job")

    def output(self):
        try:
            with open(self.file("orca.out")) as f:
                return f.readlines()
        except OSError:
            return None

    def converged(self):
        """"SCF CONVERGED" twice after a spin flip (high spin, then broken symmetry), else once."""
        lines = self.output()
        if lines is None:
            return False, "no orca.out"
        expected = 2 if "F L I P    S P I N" in "".join(lines) else 1
        found = sum("SCF CONVERGED" in line for line in lines)
        return found == expected, f"{found}/{expected} SCF"

    def backup(self):
        """For an unconverged run: delete scratch files, archive orca.inp/out/gbw as run_N.*, and
        link dft.gbw to run_N.gbw."""
        if self.converged()[0]:
            print(f"{self.config}: converged, no backup")
            return
        for fname in os.listdir(self.path):
            if ".tmp" in fname or ".bas" in fname:
                os.remove(self.file(fname))
        if not os.path.exists(self.file("orca.gbw")):
            return
        n = max(self.runs(), default=0) + 1
        subprocess.run(["cp", self.file("orca.inp"), self.file(f"run_{n}.inp")])
        subprocess.run(["cp", self.file("orca.out"), self.file(f"run_{n}.out")])
        os.rename(self.file("orca.gbw"), self.file(f"run_{n}.gbw"))
        subprocess.run(["ln", "-sf", f"run_{n}.gbw", self.file("dft.gbw")])
        print(f"{self.config}: archived orca.inp, .out and .gbw as run_{n}.*; dft.gbw -> run_{n}.gbw")

    def magnetic_moments(self):
        """Loewdin spin populations of the magnetic sites from the last such section of orca.out."""
        moments, in_section = [], False
        for line in self.output() or []:
            if "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS" in line:
                moments, in_section = [], True     # keep only the last section
                continue
            if in_section:
                if line.strip() == "" or line.startswith("----"):
                    continue
                if len(moments) < len(self.config):
                    parts = line.split()
                    moments.append({"atomic_index": int(parts[0]), "element": parts[1].rstrip(":"),
                                    "spin_population": round(float(parts[3]), 2)})
                else:
                    in_section = False
        return moments

    def check_magnetic_moments(self, tolerance=0.3):
        """Whether every site's spin population has the sign of the state ('u' positive) and the
        size of expected_moments within ``tolerance``; writes magnetic_moments.csv."""
        if not self.converged()[0]:
            print(f"{self.config}: not converged")
            return False
        moments = self.magnetic_moments()
        correct = True
        for moment in moments:
            i, population = moment["atomic_index"], moment["spin_population"]
            site = f"{self.config}: site {i} (0-based, {moment['element']})"
            if population * (1 if self.state[i] == "u" else -1) < 0:
                correct = False
                print(f"{site} has spin population {population} of the wrong sign.")
            expected = self.jobs.expected_moments[i]
            if abs(abs(population) - expected) > tolerance:
                correct = False
                print(f"{site} has spin population {population}, expected {expected} "
                      f"(tolerance: {tolerance} μB).")
        with open(self.file("magnetic_moments.csv"), "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["atomic_index", "element", "spin_population"])
            writer.writeheader()
            writer.writerows(moments)
        return correct

    def energy(self, tolerance=0.3):
        """Final single point energy (Ha) if converged with the expected moments, else None."""
        if not self.check_magnetic_moments(tolerance):
            return None
        for line in reversed(self.output()):
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1])
        return None


class CollinearOrcaJobs(CollinearJobs):
    """configs          'u'/'d' strings; the all-'u' one is the high-spin run hs_dir
    spins            spin-only |Ms| of each site, for FlipSpin/FinalMs and the multiplicity
    effective_spins  spins of the sites in energies.dat
    expected_moments expected |Loewdin spin population| of each site, checked within
                     magmom_tolerance before an energy is used
    n_cpus, mem_per_core (MB)   filled into the inputs
    """

    unit = "Ha"

    def __init__(self, configs, spins, effective_spins, expected_moments, n_cpus, mem_per_core,
                 hs_dir=None, magmom_tolerance=0.3, root=None):
        super().__init__(configs, effective_spins, root)
        self.spins = list(spins)
        self.expected_moments = list(expected_moments)
        self.n_cpus = n_cpus
        self.mem_per_core = mem_per_core
        self.hs_dir = "u" * len(self.spins) if hs_dir is None else hs_dir
        self.magmom_tolerance = magmom_tolerance

    def job(self, i, restart=False):
        return OrcaJob(self, self.configs[i], restart)

    def write_inputs(self, i):
        self.job(i).write_inputs()

    def states(self):
        return [spin_state(config) for config in self.configs]

    def energy(self, i):
        return self.job(i).energy(self.magmom_tolerance)

    @action("energies", "write energies.dat [magmom tolerance] for `spinh fit isotropic-exchange`")
    def write_energies(self, tolerance=None):
        if tolerance is not None:
            self.magmom_tolerance = float(tolerance)
        return super().write_energies()

    @action("magmoms", "check the Loewdin spin populations against expected_moments [tolerance]")
    def check_magmoms(self, tolerance=None):
        tolerance = self.magmom_tolerance if tolerance is None else float(tolerance)
        for i in range(self.n_conf):
            self.job(i).check_magnetic_moments(tolerance)

    @action("backup", "archive unconverged runs as run_N.* and link dft.gbw to their orbitals")
    def backup(self):
        for i in range(self.n_conf):
            self.job(i).backup()

    @action("restart-inputs", "rewrite the inputs for a restart from dft.gbw (no FlipSpin/FinalMs)")
    def write_restart_inputs(self):
        for i in range(self.n_conf):
            self.job(i, restart=True).write_inputs()
