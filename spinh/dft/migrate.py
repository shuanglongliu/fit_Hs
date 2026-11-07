"""``spinh migrate DIR``: the input files of a job directory from its old data.py.

The inputs used to be Python strings in data.py; they are now plain files next to main.py. This
writes them (never overwriting a file) and prints the settings of data.py that now go into
main.py, whose job set is built as in the templates (``spinh new site_vasp ...``).
"""

import os
import runpy

FILES = {
    "vasp": [("incar", "INCAR"), ("poscar", "POSCAR"), ("kpoints", "KPOINTS"), ("job_script", "vasp.job")],
    "orca": [("input_file_hs", "orca_hs.inp"), ("input_file_bs", "orca_bs.inp"), ("job_script", "orca.job")],
}

# data.py name -> where it goes now
SETTINGS = {
    "wave_dir": "SaxisJobs(wave_dir=...)",
    "emat_file": "jobs.set_frame(...)",
    "emats_file": "jobs.set_frames({site: ...})",
    "configs": "CollinearOrcaJobs(configs=...)",
    "spin_only_Szs": "CollinearOrcaJobs(spins=...)",
    "effective_spins": "CollinearOrcaJobs(effective_spins=...)",
    "expected_magmoms": "CollinearOrcaJobs(expected_moments=...)",
    "n_cpus": "CollinearOrcaJobs(n_cpus=...)",
    "mem_per_core": "CollinearOrcaJobs(mem_per_core=...)",
    "hs_dir": "CollinearOrcaJobs(hs_dir=...), relative to the directory",
}


def migrate(directory):
    directory = os.path.abspath(directory)
    data_py = os.path.join(directory, "data.py")
    if not os.path.exists(data_py):
        raise FileNotFoundError(f"No data.py in {directory}.")
    data = runpy.run_path(data_py)
    code = "orca" if "input_file_hs" in data else "vasp"
    for attr, name in FILES[code]:
        if attr not in data:
            print(f"data.py has no {attr}; {name} not written")
            continue
        path = os.path.join(directory, name)
        if os.path.exists(path):
            print(f"{name} exists, not overwritten")
            continue
        with open(path, "w") as f:
            f.write(data[attr])
        print(f"wrote {name} from data.{attr}")
    for attr, target in SETTINGS.items():
        if attr in data:
            value = data[attr]
            if attr == "hs_dir":
                value = os.path.relpath(value, directory)
            text = repr(value)
            print(f"\n{attr} -> {target}:\n  {text if len(text) < 300 else text[:300] + ' ...'}")
    print("\nmain.py now builds the job set as in the templates (spinh new <template> <dir> shows one); "
          "data.py is no longer read.")
