"""Collinear (broken-symmetry) VASP jobs. A configuration is a string of 'u' (up) and 'd'
(down), one letter per magnetic site, and is also the directory name. The INCAR template has the
fields {magmom}, {nupdown}, {U} and {J}."""

import os

from spinh.dft.slurm import run_action
from spinh.dft.vasp import VaspJob
from spinh.tables import write_collinear_energies


class CollinearVaspJobs:

    def __init__(self, data, configs, magmoms, effective_spins):
        self.data = data
        self.configs = list(configs)
        self.magmoms = magmoms
        self.effective_spins = effective_spins

    def job(self, config):
        return VaspJob(self.data, config)

    def signed(self, config, values):
        return [-v if c == "d" else v for c, v in zip(config, values)]

    def setup_jobs(self):
        """Create the directories and write the inputs."""
        for config in self.configs:
            moments = self.signed(config, self.magmoms)
            job = self.job(config)
            job.create_dir()
            job.write_inputs(magmom="".join(" {:4.1f}".format(m) for m in moments) + " 1000*0.0",
                             nupdown=sum(moments), U=0.0, J=0.0)

    def submit_jobs(self):
        for config in self.configs:
            self.job(config).submit()

    def check_convergences(self, restart=False):
        return [self.job(config).check_convergence(restart) for config in self.configs]

    def check_magmoms(self):
        """Print the site moments of the last "magnetization (x)" block of each OUTCAR."""
        for config in self.configs:
            moments = []
            try:
                with open(self.job(config).file("OUTCAR")) as f:
                    lines = f.readlines()
                for i, line in enumerate(lines):
                    if "magnetization (x)" in line:
                        moments = [float(lines[i + 4 + k].split()[-1]) for k in range(len(config))]
            except (OSError, IndexError, ValueError):
                pass
            print(f"{config}: {moments}")

    def get_energies(self):
        """Write energies.dat (eV) for the isotropic-exchange fit."""
        energies = [self.job(config).energy() for config in self.configs]
        write_collinear_energies(os.path.join(self.data.root_dir, "energies.dat"), self.configs,
                                 self.effective_spins, energies, "eV")
        return energies

    def run(self, argv=None):
        actions = {"setup": self.setup_jobs, "submit": self.submit_jobs, "check": self.check_convergences,
                   "resubmit": lambda: self.check_convergences(restart=True), "magmoms": self.check_magmoms,
                   "energies": self.get_energies}
        docs = {"setup": "create the directories and write the inputs", "submit": "submit the jobs",
                "check": "print the convergence of every job", "resubmit": "resubmit the unconverged jobs",
                "magmoms": "print the site moments", "energies": "write energies.dat"}
        return run_action(actions, argv, docs, default="check")
