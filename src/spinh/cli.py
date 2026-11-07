"""Command line: ``spinh <command> ...``; ``spinh -h`` lists the commands."""

import argparse
import os
import shutil
import sys

import yaml

TEMPLATES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "templates")


def templates():
    return sorted(d for d in os.listdir(TEMPLATES) if os.path.isdir(os.path.join(TEMPLATES, d)))


def cmd_new(args):
    if args.template not in templates():
        sys.exit(f"Unknown template {args.template!r}; choose from: {', '.join(templates())}")
    if os.path.exists(args.dir) and os.listdir(args.dir):
        sys.exit(f"{args.dir} exists and is not empty.")
    shutil.copytree(os.path.join(TEMPLATES, args.template), args.dir, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("README.md", "__pycache__"))
    print(f"Created {args.dir} from templates/{args.template}:")
    for name in sorted(os.listdir(args.dir)):
        print("  " + name)


def cmd_isotropic_exchange(args):
    from spinh.models import isotropic_exchange
    isotropic_exchange.fit(factor=args.factor, energy_threshold=args.threshold)


def cmd_anisotropic_exchange(args):
    from spinh.models import anisotropic_exchange as ae
    if args.task == "fit":
        ae.fit(factor=args.factor, quantum_correction=args.quantum_correction, iformat=args.format)
    elif args.task == "average":
        ae.average_halves()
    elif args.task == "entries":
        ae.energy_pair_entries(*args.sites, n_pair=args.n_pair)
    elif args.task == "frame":
        ae.pair_frame(args.structure, *args.sites)


def cmd_anisotropy(args):
    from spinh.models import single_ion_anisotropy as sia
    with open(args.config) as f:
        config = yaml.safe_load(f)
    cfg = dict(config[args.task])
    cfg.setdefault("spin", config.get("spin"))
    sia.TASKS[args.task](cfg)


def cmd_polarization(args):
    from spinh.models import polarization
    polarization.MODELS[args.model]()


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spinh", description="Spin Hamiltonian fitting.")
    sub = parser.add_subparsers(dest="command", required=True)

    def add(name, func, help):
        p = sub.add_parser(name, help=help, description=help)
        p.set_defaults(func=func)
        if name != "new":
            p.add_argument("-d", "--dir", default=".", help="working directory (default: .)")
        return p

    p = add("new", cmd_new, "copy a template into a new directory")
    p.add_argument("template", help="one of: " + ", ".join(templates()))
    p.add_argument("dir")

    p = add("isotropic-exchange", cmd_isotropic_exchange,
            "fit isotropic J from energies.dat and pairs.yaml")
    p.add_argument("--factor", type=float, default=-2.0, help="H = factor sum J S_i.S_j (default -2)")
    p.add_argument("--threshold", type=float, default=None,
                   help="drop configurations more than this many cm^-1 above the lowest")

    p = add("anisotropic-exchange", cmd_anisotropic_exchange, "3x3 exchange matrices (input.yaml)")
    p.add_argument("task", choices=["fit", "average", "entries", "frame"],
                   help="fit: fit J; average: energies.dat -> energies_average.dat (row i with "
                        "row i + n/2); entries: print input.yaml entries for two sites; "
                        "frame: pair reference frame from a structure file")
    p.add_argument("sites", nargs="*", type=int, help="two site indices (entries, frame)")
    p.add_argument("--factor", type=float, default=-2.0, help="H = factor sum S_i^T J S_j (default -2)")
    p.add_argument("--quantum-correction", action="store_true", help="Ruiz et al. (2003) correction")
    p.add_argument("--format", type=int, default=1, choices=[1, 2], help="1: per pair; 2: one line")
    p.add_argument("--n-pair", type=int, default=9, help="entries: row k pairs with row k + n_pair")
    p.add_argument("--structure", default="POSCAR", help="frame: structure file read by ase")

    p = add("anisotropy", cmd_anisotropy, "single-ion anisotropy tasks configured in anisotropy.yaml")
    p.add_argument("task", choices=["average", "change-frame", "mae", "plot", "fit-DE", "fit-D-tensor",
                                    "fit-curve", "fit-surface"])
    p.add_argument("--config", default="anisotropy.yaml")

    p = add("polarization", cmd_polarization, "fit the spin-driven polarization (input.yaml, dipoles.dat)")
    p.add_argument("model", choices=["alpha-pi", "full-p-v1", "full-p-v2"])

    args = parser.parse_args(argv)
    if args.command == "anisotropic-exchange" and args.task in ("entries", "frame") and len(args.sites) != 2:
        parser.error(f"{args.task} needs two site indices")
    if args.command != "new":
        os.chdir(args.dir)
    args.func(args)


if __name__ == "__main__":
    main()
