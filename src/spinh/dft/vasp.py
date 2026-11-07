"""One VASP run in one directory.

``data`` is the job directory's data.py module: root_dir, incar (a format string), poscar,
kpoints and job_script. POTCAR is linked from root_dir.
"""

import os
import subprocess

import numpy as np

from spinh.dft.slurm import sbatch


class VaspJob:

    def __init__(self, data, dir_name):
        self.data = data
        self.dir_name = dir_name
        self.path = os.path.join(data.root_dir, dir_name)

    def file(self, name):
        return os.path.join(self.path, name)

    def create_dir(self):
        os.makedirs(self.path, exist_ok=True)

    def write_inputs(self, wave_dir=None, **fields):
        """INCAR (``incar.format(**fields)``), POSCAR, KPOINTS, vasp.job and a POTCAR link; a
        WAVECAR link to ``wave_dir`` unless this is that directory."""
        print(self.path)
        for name, text in (("INCAR", self.data.incar.format(**fields)), ("POSCAR", self.data.poscar),
                           ("KPOINTS", self.data.kpoints), ("vasp.job", self.data.job_script)):
            with open(self.file(name), "w") as f:
                f.write(text)
        self.link(os.path.join(self.data.root_dir, "POTCAR"), "POTCAR")
        if wave_dir is not None and self.dir_name != wave_dir:
            self.link(os.path.join(self.data.root_dir, wave_dir, "WAVECAR"), "WAVECAR")

    def link(self, target, name):
        subprocess.run(["ln", "-sf", target, self.file(name)])

    def submit(self):
        sbatch(self.path, "vasp.job")

    def converged(self, criterion="ediff", de0=1.e-8):
        """(converged, last energy change as a string).

        criterion "ediff": OUTCAR says "EDIFF is reached"; the change is read from OSZICAR.
        criterion "output": the last line of the job's ``output`` with a ':' has |dE| < de0.
        """
        try:
            with open(self.file("OSZICAR" if criterion == "ediff" else "output")) as f:
                de = [line for line in f if ":" in line][-1].split()[3]
        except (OSError, IndexError):
            return False, "NA"
        if criterion == "ediff":
            try:
                with open(self.file("OUTCAR")) as f:
                    return "EDIFF is reached" in f.read(), de
            except OSError:
                return False, de
        try:
            return abs(float(de)) < de0, de
        except ValueError:
            return False, de

    def check_convergence(self, restart=False, criterion="ediff", de0=1.e-8):
        """Print the convergence; with ``restart``, resubmit if not converged."""
        converged, de = self.converged(criterion, de0)
        print("dir_name = {:>30s} , de = {:>15s} , convergence = {:d}.".format(self.dir_name, de, converged))
        if restart and not converged:
            self.submit()
        return converged

    def energy(self):
        """Last energy(sigma->0) in OUTCAR (eV), 0 if there is none."""
        energy = 0
        try:
            with open(self.file("OUTCAR")) as f:
                for line in f:
                    if "energy(sigma->0)" in line:
                        energy = float(line.split()[-1])
        except (OSError, ValueError):
            pass
        return energy

    def dipole(self, crystal=False):
        """[px, py, pz, |p|] from OUTCAR: the dipole correction ("dipolmoment") or, with
        ``crystal``, the Berry-phase electronic dipole."""
        tag, m, n = ("electronic dipole moment", 5, 8) if crystal else ("dipolmoment", 1, 4)
        try:
            with open(self.file("OUTCAR")) as f:
                for line in f:
                    if tag in line:
                        p = [float(v) for v in line.split()[m:n]]
                        return p + [np.linalg.norm(p)]
        except (OSError, ValueError):
            pass
        return [np.nan] * 4

    def local_moments(self, n_site, kind="magnetization"):
        """(n_site, 3) site moments from the first complete x, y, z block of OUTCAR; kind is
        "magnetization" or "orbital moment". NaN if missing."""
        try:
            with open(self.file("OUTCAR")) as f:
                lines = f.readlines()
            start = {}
            for i, line in enumerate(lines):
                for axis in "xyz":
                    if f"{kind} ({axis})" in line:
                        start[axis] = i
                if "z" in start:
                    break
            return np.array([[float(lines[start[axis] + 4 + i].split()[-1]) for axis in "xyz"]
                             for i in range(n_site)])
        except (OSError, KeyError, IndexError, ValueError):
            return np.full((n_site, 3), np.nan)

    def occupation_matrices(self):
        """Write occmat (occupation matrices) and occmat.eigen (their eigenvalues and
        eigenvectors) from the last iteration in OUTCAR of a noncollinear run."""
        with open(self.file("OUTCAR")) as f:
            lines = f.readlines()
        start = max([i for i, line in enumerate(lines) if "Iteration" in line], default=0)
        lines = lines[start:]
        atoms = [(i, int(line.split()[-1])) for i, line in enumerate(lines) if "atom =" in line and "l =" in line]
        n_spin = 4

        entry = "# {:3d}\n".format(len(atoms))
        eigen = ""
        for i_atom, (i, l) in enumerate(atoms):
            tag = "# {:3d} {:3d} {:3d}\n".format(i_atom + 1, l, n_spin)
            entry += tag
            entry += "".join(line for line in lines[i + 4:i + 4 + n_spin * (2 * l + 1 + 3)] if line.strip())
            entry += "# \n"
            eigen += tag
            for line in lines[i + 38:i + 38 + 2 * (2 * l + 1)]:
                if line.strip():
                    o, v = line.replace("o =", "").strip().split("v =")
                    eigen += v + "   #   " + o + "\n"
        with open(self.file("occmat"), "w") as f:
            f.write(entry.replace("spin", "# spin"))
        with open(self.file("occmat.eigen"), "w") as f:
            f.write(eigen)
        print(f"{self.path}: wrote occmat and occmat.eigen")
