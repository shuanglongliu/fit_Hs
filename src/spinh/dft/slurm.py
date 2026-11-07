"""Slurm submission and command-line actions for job sets."""

import subprocess
import sys
import time


def sbatch(path, script):
    print(path)
    subprocess.run(["sbatch", script], cwd=path)
    time.sleep(1)


def run_action(actions, argv=None, docs=None, default="print"):
    """Call ``actions[name](*args)`` for ``argv = [name, *args]`` (default: sys.argv[1:])."""
    argv = sys.argv[1:] if argv is None else argv
    name, args = (argv[0], argv[1:]) if argv else (default, [])
    if name not in actions:
        docs = docs or {}
        width = max(len(n) for n in actions)
        lines = "\n".join(f"  {n:{width}s}  {docs.get(n, '')}" for n in actions)
        sys.exit(f"usage: python main.py ACTION [ARGS]   (default: {default})\n\nactions:\n{lines}"
                 + ("" if name in ("-h", "--help", "help") else f"\n\nunknown action: {name}"))
    return actions[name](*args)
