"""Plots shared by the fits."""

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
mpl.rcParams["font.size"] = 12


def plot_fit(y, y_fit, fname, xlabel="", ylabel="", title="", x=None, x_fit=None):
    """DFT values (dots) and fitted values against x: squares, or a line if ``x_fit`` is given."""
    if x is None:
        x = range(1, len(y) + 1)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.set_title(title, {"fontsize": 10})
    ax.plot(x, y, "bo", label="DFT")
    if x_fit is None:
        ax.plot(x, y_fit, "rs", label="Fit")
    else:
        ax.plot(x_fit, y_fit, "r-", label="Fit")
    ax.legend(frameon=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close(fig)
