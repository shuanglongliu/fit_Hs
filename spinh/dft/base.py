"""What every set of DFT jobs has in common: one directory per configuration under a root
directory, and the actions run by ``python main.py ACTION [ARGS]``.

The root directory is the directory of main.py unless given: it holds the input files shared by
the jobs (INCAR, POSCAR, ... or the ORCA inputs) and the job directories.

A subclass provides
    dir_names              the job directories, one per configuration
    job(i)                 the job of configuration i (VaspJob or OrcaJob)
    write_inputs(i)        write the input files of configuration i
    print_configurations   the "print" action
    write_energies         the "energies" action
and adds its own actions with the @action decorator.
"""

import os
import sys


def script_dir():
    """Directory of the script being run (main.py), or the working directory."""
    path = getattr(sys.modules.get("__main__"), "__file__", None)
    return os.path.dirname(os.path.abspath(path)) if path else os.getcwd()


def action(name, help):
    """Mark a method as the action ``python main.py <name> [ARGS]``; ``help`` is listed by
    ``python main.py help``. ARGS are passed to the method as strings."""
    def mark(method):
        method.action = (name, help)
        return method
    return mark


class JobSet:

    default_action = "print"

    def __init__(self, root=None):
        self.root = script_dir() if root is None else os.path.abspath(root)

    def path(self, name):
        """A file or directory under root."""
        return os.path.join(self.root, name)

    @property
    def n_conf(self):
        return len(self.dir_names)

    def job(self, i):
        raise NotImplementedError

    def write_inputs(self, i):
        raise NotImplementedError

    def converged(self, i):
        """(converged, information such as the last energy change) of job i."""
        return self.job(i).converged()

    # -- actions ----------------------------------------------------------------------------------

    @action("print", "print the configurations and their directories")
    def print_configurations(self):
        raise NotImplementedError

    @action("setup", "create the directories and write the inputs")
    def setup(self):
        for i in range(self.n_conf):
            self.job(i).create_dir()
            self.write_inputs(i)

    @action("submit", "submit every job")
    def submit(self):
        for i in range(self.n_conf):
            self.job(i).submit()

    @action("check", "print the convergence of every job")
    def check(self):
        result = []
        for i in range(self.n_conf):
            converged, info = self.converged(i)
            print("{:>30s}  {:>15s}  {}".format(self.dir_names[i], info, "converged" if converged else "NOT converged"))
            result.append(converged)
        return result

    @action("resubmit", "resubmit the jobs that are not converged")
    def resubmit(self):
        for i, converged in enumerate(self.check()):
            if not converged:
                self.job(i).submit()

    @action("energies", "write energies.dat")
    def write_energies(self):
        raise NotImplementedError

    # -- command line -----------------------------------------------------------------------------

    def actions(self):
        """{name: (bound method, help)}, base-class actions first."""
        found = {}
        for cls in reversed(type(self).__mro__):
            for attr, member in vars(cls).items():
                if hasattr(member, "action"):
                    name, help = member.action
                    found[name] = (getattr(self, attr), help)
        return found

    def usage(self):
        actions = self.actions()
        width = max(len(name) for name in actions)
        lines = "\n".join(f"  {name:{width}s}  {help}" for name, (_, help) in actions.items())
        return f"usage: python main.py ACTION [ARGS]   (default: {self.default_action})\n\nactions:\n{lines}"

    def run(self, argv=None):
        """Run the action named in ``argv`` (default: sys.argv[1:])."""
        argv = sys.argv[1:] if argv is None else list(argv)
        name, args = (argv[0], argv[1:]) if argv else (self.default_action, [])
        actions = self.actions()
        if name in ("-h", "--help", "help"):
            print(self.usage())
            return None
        if name not in actions:
            sys.exit(self.usage() + f"\n\nunknown action: {name}")
        return actions[name][0](*args)
