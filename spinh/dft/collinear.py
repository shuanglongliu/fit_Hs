"""Collinear (broken-symmetry) jobs. A configuration is a string of 'u' (up) and 'd' (down), one
letter per magnetic site, and is also the directory name.

CollinearJobs        what VASP and ORCA share: printing the configurations and writing
                     energies.dat (signed effective spins, then the energy) for
                     ``spinh fit isotropic-exchange``
CollinearVaspJobs    VASP; the INCAR template has the fields {magmom}, {nupdown}, {U} and {J}
CollinearOrcaJobs    ORCA, in spinh.dft.orca
"""

from spinh.dft.base import JobSet, action
from spinh.dft.vasp import VaspJob
from spinh.tables import write_collinear_energies


def signed(config, values):
    """values with the sign of each site's spin: 'u' +, 'd' -."""
    return [-v if c == "d" else v for c, v in zip(config, values)]


class CollinearJobs(JobSet):

    unit = "eV"     # of the energies

    def __init__(self, configs, effective_spins, root=None):
        super().__init__(root)
        self.configs = list(configs)
        self.effective_spins = list(effective_spins)
        for config in self.configs:
            if len(config) != len(self.effective_spins) or set(config) - set("ud"):
                raise ValueError(f"Configuration {config!r}: expected {len(self.effective_spins)} letters 'u' or 'd'.")

    @property
    def dir_names(self):
        return self.configs

    def states(self):
        """The spin state of each configuration written to energies.dat."""
        return self.configs

    def energy(self, i):
        """Energy of configuration i, None if there is none."""
        raise NotImplementedError

    @action("print", "print the configurations and the signed effective spins")
    def print_configurations(self):
        for config in self.configs:
            print(f"{config:20s}" + "".join(f" {s:5.1f}" for s in signed(config, self.effective_spins)))

    @action("energies", "write energies.dat for `spinh fit isotropic-exchange`")
    def write_energies(self):
        energies = [self.energy(i) for i in range(self.n_conf)]
        for config, energy in zip(self.configs, energies):
            if energy is None:
                print(f"{config}: no energy, left out")
        fname = self.path("energies.dat")
        write_collinear_energies(fname, self.states(), self.effective_spins, energies, self.unit)
        print(f"Wrote {fname} ({self.unit}).")
        return energies


class CollinearVaspJobs(CollinearJobs):

    def __init__(self, configs, magmoms, effective_spins, root=None):
        super().__init__(configs, effective_spins, root)
        self.magmoms = list(magmoms)

    def job(self, i):
        return VaspJob(self.root, self.configs[i])

    def write_inputs(self, i):
        moments = signed(self.configs[i], self.magmoms)
        self.job(i).write_inputs(magmom="".join(" {:4.1f}".format(m) for m in moments) + " 1000*0.0",
                                 nupdown=sum(moments), U=0.0, J=0.0)

    def energy(self, i):
        energy = self.job(i).energy()
        return None if energy == 0 else energy

    @action("magmoms", "print the site moments of every job")
    def print_magmoms(self):
        """The site moments of the last "magnetization (x)" block of each OUTCAR."""
        for i, config in enumerate(self.configs):
            moments = []
            try:
                with open(self.job(i).file("OUTCAR")) as f:
                    lines = f.readlines()
                for k, line in enumerate(lines):
                    if "magnetization (x)" in line:
                        moments = [float(lines[k + 4 + s].split()[-1]) for s in range(len(config))]
            except (OSError, IndexError, ValueError):
                pass
            print(f"{config}: {moments}")
