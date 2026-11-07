"""Regression tests: python -m unittest discover tests

The files in tests/reference were produced by the scripts this package replaced.
"""

import contextlib
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import yaml

from spinh import cli, geometry
from spinh.models import single_ion_anisotropy
from spinh.sphere import Sphere
from spinh.spin_operators import spin_operators

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(HERE, "reference")


def reference(name):
    with open(os.path.join(REFERENCE, name)) as f:
        return f.read()


def snapshot(root):
    """Files (hash) and links (target) under root, like the stored reference trees."""
    root = os.path.abspath(root)
    out = []
    for d, dirs, files in sorted(os.walk(root)):
        dirs.sort()
        for f in sorted(files + [x for x in dirs if os.path.islink(os.path.join(d, x))]):
            p = os.path.join(d, f)
            rel = os.path.relpath(p, root)
            if rel.endswith(".py") or "__pycache__" in rel:
                continue
            if os.path.islink(p):
                out.append(f"L {rel} -> {os.readlink(p).replace(root, '<ROOT>')}")
            else:
                with open(p, "rb") as fh:
                    out.append(f"F {rel} {hashlib.sha256(fh.read()).hexdigest()[:16]}")
    return out


class TempDirTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cwd = os.getcwd()

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.tmp)

    def new(self, template):
        path = os.path.join(self.tmp, template)
        with contextlib.redirect_stdout(io.StringIO()):
            cli.main(["new", template, path])
        return path

    def spinh(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            cli.main(list(argv))
        os.chdir(self.cwd)
        return out.getvalue()


class TestGeometry(unittest.TestCase):

    def test_rotate_frame(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            alpha, beta, gamma = rng.uniform(0, 2 * np.pi, 3)
            emat = geometry.rotate_frame(np.eye(3), alpha, beta, gamma)
            # the new z axis points along (beta, gamma)
            np.testing.assert_allclose(emat[2], geometry.sph2cart([1, beta, gamma], rad=True), atol=1e-12)
            # alpha = 0: tilt by beta, turn x and y about the new z by gamma (old fit_D_and_E: y then new z)
            # gamma = 0: as above with alpha; the same frame as rotating about the old z after the tilt
            Ry = geometry.euler_zyz(0, beta, 0, degrees=False)
            Rz = geometry.euler_zyz(alpha, 0, 0, degrees=False)
            np.testing.assert_allclose(geometry.rotate_frame(np.eye(3), alpha, beta, 0), (Ry @ Rz).T, atol=1e-12)
            Rz = geometry.euler_zyz(gamma, 0, 0, degrees=False)
            np.testing.assert_allclose(geometry.rotate_frame(np.eye(3), 0, beta, gamma), (Rz @ Ry).T, atol=1e-12)

    def test_sphere(self):
        sphere = Sphere()
        for index in [0, 1, 3599, 3600, 3601, 4319101, sphere.n_direction - 1]:
            self.assertEqual(sphere.index(sphere.direction(index)), index)

    def test_spin_operators(self):
        for S in (0.5, 1, 1.5, 4):
            o = spin_operators(S)
            np.testing.assert_allclose(o.Sx @ o.Sy - o.Sy @ o.Sx, 1j * o.Sz, atol=1e-12)
            np.testing.assert_allclose(o.Sx @ o.Sx + o.Sy @ o.Sy + o.Sz @ o.Sz, o.ss * o.ID, atol=1e-12)


class TestTemplates(TempDirTest):
    """main.py setup writes the same directories and inputs as the old job directories."""

    def run_template(self, template, *actions, edit=None):
        path = self.new(template)
        shutil.copy(os.path.join(REFERENCE, "POTCAR"), path)
        if edit:
            with open(os.path.join(path, "main.py")) as f:
                text = f.read()
            with open(os.path.join(path, "main.py"), "w") as f:
                f.write(edit(text))
        for action in actions:
            result = subprocess.run([sys.executable, "main.py", action], cwd=path, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        return path

    def check_tree(self, template, name, *actions):
        path = self.run_template(template, *actions)
        expected = json.loads(reference("job_trees.json"))[name]
        self.assertEqual(snapshot(path), expected)

    def test_site(self):
        self.check_tree("site_vasp", "site_magnetic_moment", "setup")

    def test_saxis(self):
        self.check_tree("saxis_vasp", "saxis", "setup")

    def test_collinear_vasp(self):
        self.check_tree("collinear_vasp", "collinear_vasp", "setup")

    def test_collinear_orca(self):
        self.check_tree("collinear_orca", "collinear_orca", "setup")

    def test_site_spin_directions(self):
        """The spin set: MAGMOM = moments opposite to the listed spin directions."""
        path = self.run_template("site_vasp", "setup", edit=lambda t: t.replace(
            "    myjob = do(direction_emin = (0.0, 0.0))\n    # myjob = do_spin(",
            "    # myjob = do(direction_emin = (0.0, 0.0))\n    myjob = do_spin("))
        with open(os.path.join(path, "0_0_0", "INCAR")) as f:
            incar = f.read()
        self.assertIn("MAGMOM =    0.0000000000   0.0000000000  -2.0000000000", incar)
        self.assertEqual(len([d for d in os.listdir(path) if d[0].isdigit()]), 30)

    def test_harvest(self):
        """energies.dat from the jobs' OUTCAR files."""
        path = self.run_template("site_vasp", "setup")
        names = sorted(d for d in os.listdir(path) if d[0].isdigit())
        for k, name in enumerate(names):
            with open(os.path.join(path, name, "OUTCAR"), "w") as f:
                f.write(f"  energy  without entropy=  1.0  energy(sigma->0) =  {-100 - k:.6f}\n")
        subprocess.run([sys.executable, "main.py", "all-energies"], cwd=path, capture_output=True, check=True)
        data = np.loadtxt(os.path.join(path, "energies.dat"))
        self.assertEqual(data.shape, (len(names), 21))
        self.assertEqual(sorted(data[:, -1]), sorted(-100.0 - np.arange(len(names))))


class TestFits(TempDirTest):

    def test_isotropic_exchange(self):
        path = self.new("fit_isotropic_exchange")
        self.spinh("isotropic-exchange", "-d", path)
        with open(os.path.join(path, "Js.csv")) as f:
            self.assertEqual(f.read(), reference("Js.csv"))

    def test_anisotropic_exchange(self):
        path = self.new("fit_anisotropic_exchange")
        out = self.spinh("anisotropic-exchange", "fit", "-d", path, "--format", "2")
        meV = [float(v) for v in out.split("J (meV):")[1].split("\n")[0].split()]
        self.assertEqual(len(meV), 27)
        self.assertAlmostEqual(meV[0], -1.884, places=3)
        self.assertAlmostEqual(meV[-1], -1.865, places=3)

    def test_polarization(self):
        for template, model, files in [
            ("fit_polarization_alpha_pi", "alpha-pi", {"log": "alpha_pi.log"}),
            ("fit_polarization_full_p_v1", "full-p-v1", {"P_total.log": "full_p_v1_P_total.log"}),
            ("fit_polarization_full_p_v2", "full-p-v2", {"P_for_latex.log": "full_p_v2_P_for_latex.log",
                                                         "relative_contribution.log": "full_p_v2_relative_contribution.log"})]:
            path = self.new(template)
            self.spinh("polarization", model, "-d", path)
            for out, ref in files.items():
                with open(os.path.join(path, out)) as f:
                    self.assertEqual(f.read(), reference(ref), out)

    def test_anisotropy_tasks(self):
        path = self.new("fit_anisotropy")
        for task in ("average", "change-frame", "mae", "plot", "fit-DE", "fit-D-tensor", "fit-curve", "fit-surface"):
            self.spinh("anisotropy", task, "-d", path)
        for out, ref in (("energies_average.dat", "energies_average.dat"), ("new_frame.dat", "new_frame.dat")):
            with open(os.path.join(path, out)) as f:
                self.assertEqual(f.read(), reference(ref), out)

    def test_fit_DE_general_frame(self):
        """All three Euler angles free: D, E and the eigenframe are recovered."""
        D, E, angles = -0.8, -0.15, np.radians([25, 30, 60])
        frame = geometry.rotate_frame(np.eye(3), *angles)
        os.chdir(self.tmp)
        with open("e.dat", "w") as f:
            for th in range(0, 181, 15):
                for ph in ([0] if th in (0, 180) else range(0, 360, 20)):
                    n = frame @ geometry.sph2cart([1, th, ph])
                    e = 4 * (D * n[2]**2 + E * (n[0]**2 - n[1]**2)) / 1000 - 50
                    f.write(f"{th} {ph} {e:.12f}\n")
        cfg = {"input": "e", "columns": [1, 2, 3], "spin": 2, "model": "classical",
               "bounds": {"D": [-2, 0], "E": [-0.5, 0], "E0": [-10, 10], "alpha": [0, 90], "beta": [0, 90],
                          "gamma": [0, 180]},
               "p0": {"D": -0.5, "E": -0.1, "E0": 0, "alpha": 20, "beta": 20, "gamma": 40}}
        with contextlib.redirect_stdout(io.StringIO()) as out:
            single_ion_anisotropy.fit_DE(cfg)
        values = out.getvalue().split("\n")[0].replace(",", "").split()
        self.assertAlmostEqual(float(values[2]), D, places=4)
        self.assertAlmostEqual(float(values[6]), E, places=4)
        fitted = np.loadtxt("d_spin-orbit.dat")[:, 1:]
        np.testing.assert_allclose(np.abs(np.sum(fitted * frame, axis=1)), 1, atol=1e-4)


if __name__ == "__main__":
    unittest.main()
