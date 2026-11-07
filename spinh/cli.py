"""Command line: ``spinh <command> ...``; ``spinh -h`` lists the commands.

    spinh new [TEMPLATE DIR]          copy a template (no arguments: list them)
    spinh fit MODEL                   fit a model to the data named in fit.yaml
    spinh tool TOOL [ARGS]            prepare or inspect data
    spinh migrate DIR                 input files of an old job directory from its data.py

fit.yaml: the top-level keys hold the settings shared by every model and tool; a section named
after a model or tool holds its own settings and overrides them.
"""

import argparse
import os
import shutil
import sys

import yaml

TEMPLATES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

TEMPLATE_HELP = {
    "collinear_vasp": "VASP jobs for 'u'/'d' configurations (isotropic exchange)",
    "collinear_orca": "broken-symmetry ORCA jobs for 'u'/'d' configurations (isotropic exchange)",
    "site_vasp": "noncollinear VASP jobs, one moment direction per site (anisotropy, exchange matrices, dipoles)",
    "saxis_vasp": "noncollinear VASP jobs, one spin axis for the whole system (anisotropy)",
    "fit_isotropic_exchange": "spinh fit isotropic-exchange, with example data",
    "fit_anisotropic_exchange": "spinh fit anisotropic-exchange, with example data",
    "fit_anisotropy": "spinh fit anisotropy-* and the anisotropy tools, with synthetic data",
    "fit_polarization_alpha_pi": "spinh fit polarization-alpha-pi, with example data",
    "fit_polarization_full_p_v1": "spinh fit polarization-full-p-v1, with example data",
    "fit_polarization_full_p_v2": "spinh fit polarization-full-p-v2, with example data",
}


def fit_models():
    """{model: (function of the settings, help)}."""
    from spinh.models import anisotropic_exchange, isotropic_exchange, polarization
    from spinh.models import single_ion_anisotropy as sia
    return {
        "isotropic-exchange": (isotropic_exchange.fit, "isotropic J from collinear energies"),
        "anisotropic-exchange": (anisotropic_exchange.fit, "3x3 exchange matrices from noncollinear energy pairs"),
        "anisotropy-DE": (sia.fit_DE, "single-ion D, E and the eigenframe"),
        "anisotropy-D-tensor": (sia.fit_D_tensor, "single-ion symmetric D tensor"),
        "anisotropy-curve": (sia.fit_curve, "a (theta - b)^2 + c along one cut"),
        "anisotropy-surface": (sia.fit_surface, "a theta'^2 + c about a fitted axis"),
        "polarization-alpha-pi": (polarization.fit_alpha_pi, "p from alpha (cross) and pi (dot) terms"),
        "polarization-full-p-v1": (polarization.fit_full_p_v1, "full P tensors, symmetry as constraints"),
        "polarization-full-p-v2": (polarization.fit_full_p_v2, "full P tensors of the irreducible pairs"),
    }


def config_tools():
    """Tools configured in fit.yaml: {tool: (function of the settings, help)}."""
    from spinh.models import single_ion_anisotropy as sia
    return {
        "average-opposite": (sia.average, "average E(n) and E(-n) -> <data>_average.dat"),
        "change-frame": (sia.change_frame, "directions in a frame with z along z_direction -> new_frame.dat"),
        "mae": (sia.mae, "out-of-plane and in-plane anisotropy energies"),
        "plot": (sia.polar_plot, "polar map of the energies -> <data>.pdf"),
    }


def templates():
    return sorted(d for d in os.listdir(TEMPLATES) if os.path.isdir(os.path.join(TEMPLATES, d))
                  and not d.startswith("__"))


def load_config(path, name):
    """The top-level settings of ``path`` (other than the sections of models and tools),
    updated with the section ``name``."""
    if not os.path.exists(path):
        sys.exit(f"{os.path.abspath(path)} not found (spinh new shows the fit templates).")
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    sections = set(fit_models()) | set(config_tools())
    cfg = {k: v for k, v in data.items() if k not in sections}
    cfg.update(data.get(name) or {})
    return cfg


def cmd_new(args):
    if args.template is None:
        width = max(len(t) for t in templates())
        print("spinh new TEMPLATE DIR, with TEMPLATE one of:\n")
        for t in templates():
            print(f"  {t:{width}s}  {TEMPLATE_HELP.get(t, '')}")
        return
    if args.template not in templates():
        sys.exit(f"Unknown template {args.template!r}; choose from: {', '.join(templates())}")
    if args.dir is None:
        sys.exit("spinh new TEMPLATE DIR: give the directory to create.")
    if os.path.exists(args.dir) and os.listdir(args.dir):
        sys.exit(f"{args.dir} exists and is not empty.")
    shutil.copytree(os.path.join(TEMPLATES, args.template), args.dir, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("README.md", "__pycache__"))
    print(f"Created {args.dir} from the template {args.template}:")
    for name in sorted(os.listdir(args.dir)):
        print("  " + name)


def cmd_fit(args):
    fit_models()[args.model][0](load_config(args.config, args.model))


def cmd_tool(args):
    from spinh.models import anisotropic_exchange as ae
    if args.tool in config_tools():
        config_tools()[args.tool][0](load_config(args.config, args.tool))
    elif args.tool == "average-halves":
        ae.average_halves(args.input, args.output)
    elif args.tool == "entries":
        ae.energy_pair_entries(args.i, args.j, args.input, n_pair=args.n_pair)
    elif args.tool == "pair-frame":
        ae.pair_frame(args.structure, args.i, args.j)


def cmd_migrate(args):
    from spinh.dft.migrate import migrate
    migrate(args.dir)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spinh", description="Spin Hamiltonian parameters from DFT.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="spinh <command> -h describes a command.")
    sub = parser.add_subparsers(dest="command", required=True, metavar="command")

    def workdir(p):
        p.add_argument("-d", "--dir", default=".", help="working directory (default: .)")

    p = sub.add_parser("new", help="copy a template into a new directory (no arguments: list the templates)")
    p.add_argument("template", nargs="?")
    p.add_argument("dir", nargs="?")
    p.set_defaults(func=cmd_new)

    models = fit_models()
    p = sub.add_parser("fit", help="fit a model; settings and data file names in fit.yaml",
                       formatter_class=argparse.RawDescriptionHelpFormatter,
                       description="Fit a model. The settings are the top-level keys of fit.yaml and its "
                                   "section named after the model.\n\nmodels:\n"
                                   + "\n".join(f"  {m:24s} {h}" for m, (_, h) in models.items()))
    p.add_argument("model", choices=list(models), metavar="model")
    p.add_argument("-c", "--config", default="fit.yaml", help="settings file (default: fit.yaml)")
    workdir(p)
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("tool", help="prepare or inspect data")
    tools = p.add_subparsers(dest="tool", required=True, metavar="tool")
    for name, (_, help) in config_tools().items():
        t = tools.add_parser(name, help=help + " (settings: fit.yaml)", description=help)
        t.add_argument("-c", "--config", default="fit.yaml", help="settings file (default: fit.yaml)")
        workdir(t)
    t = tools.add_parser("average-halves", help="energy of row k of energies.dat averaged with row k + n/2",
                         description="Average the energy of row k with row k + n/2 (the same configurations "
                                     "with the other spins reversed).")
    t.add_argument("--input", default="energies.dat")
    t.add_argument("--output", default="energies_average.dat")
    workdir(t)
    t = tools.add_parser("entries", help="fit.yaml energies entries for sites I, J: row k with row k + n_pair")
    t.add_argument("i", type=int)
    t.add_argument("j", type=int)
    t.add_argument("--n-pair", type=int, default=9)
    t.add_argument("--input", default="energies.dat")
    workdir(t)
    t = tools.add_parser("pair-frame", help="reference frame of the pair I, J from a structure (needs ase)")
    t.add_argument("i", type=int)
    t.add_argument("j", type=int)
    t.add_argument("--structure", default="POSCAR", help="structure file read by ase (default: POSCAR)")
    workdir(t)
    p.set_defaults(func=cmd_tool)

    p = sub.add_parser("migrate", help="input files of an old job directory from its data.py")
    p.add_argument("dir", nargs="?", default=".")
    p.set_defaults(func=cmd_migrate)

    args = parser.parse_args(argv)
    if hasattr(args, "dir") and args.command in ("fit", "tool"):
        os.chdir(args.dir)
    args.func(args)


if __name__ == "__main__":
    main()
