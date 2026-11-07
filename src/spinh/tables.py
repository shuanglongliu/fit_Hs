"""Writers for the data files that connect the DFT jobs and the fits."""


def write_direction_table(fname, configurations, configurations_local, rows, row_format):
    """Write one line per configuration: theta, phi of every site, then the values in ``rows``.

    The block in the global frame is commented out, so ``np.loadtxt`` reads the local block.
    """
    def lines(configurations, prefix):
        out = ""
        for conf, row in zip(configurations, rows):
            out += prefix + "".join("{:6.1f}  {:6.1f}  ".format(*d) for d in conf) + row_format.format(*row) + "\n"
        return out

    with open(fname, "w") as f:
        f.write("# global reference frame\n")
        f.write(lines(configurations, "#"))
        f.write("\n\n")
        f.write("# local reference frame\n")
        f.write(lines(configurations_local, ""))


def write_collinear_energies(fname, configs, effective_spins, energies, unit):
    """energies.dat for the isotropic-exchange fit: signed effective spins, then the energy.

    ``configs`` are 'u'/'d' strings; configurations whose energy is None are skipped.
    """
    with open(fname, "w") as f:
        f.write("# " + "   ".join(f"S{i + 1}" for i in range(len(effective_spins))) + f"   Energy ({unit})\n")
        for config, energy in zip(configs, energies):
            if energy is None:
                continue
            spins = [s if c == "u" else -s for c, s in zip(config, effective_spins)]
            f.write("".join(f" {s:5.1f}" for s in spins) + f" {energy:18.10f}\n")
