"""Linear least squares with optional linear equality constraints."""

import numpy as np
from scipy.linalg import null_space


def solve(A, y, C=None, z=None, require_determined=False):
    """Minimize |A x - y| subject to C x = z.

    x is written as x_p + N w, where x_p solves C x = z and the columns of N span the null space
    of C; w is found by least squares. When the data do not determine every free parameter, the
    minimum-norm solution is returned with a warning, or ValueError is raised if
    ``require_determined``.
    """
    A = np.asarray(A, dtype=float)
    if C is None or len(C) == 0:
        x_p = np.zeros(A.shape[1])
        N = np.eye(A.shape[1])
    else:
        x_p = np.linalg.lstsq(C, z, rcond=None)[0]
        N = null_space(C)
    n_free = N.shape[1]
    if require_determined and A.shape[0] < n_free:
        raise ValueError(f"Underdetermined system: {A.shape[0]} equation(s) for {n_free} unknown(s).")
    w, _, rank, _ = np.linalg.lstsq(A @ N, y - A @ x_p, rcond=None)
    if rank < n_free:
        message = (f"{n_free - rank} combination(s) of the {n_free} free parameters "
                   f"are not determined by the data")
        if require_determined:
            raise ValueError(message + ".")
        print(f"Warning: {message}; returning the minimum-norm solution.")
    return x_p + N @ w


def report(A, x, y, unit):
    """Print the size of the system and, when overdetermined, the residual."""
    n_eq, n_x = A.shape
    if n_eq == n_x:
        print(f"Square system: {n_eq} equations for {n_x} unknowns.")
    else:
        print(f"{n_eq} equations for {n_x} unknowns: residual |A x - y| = "
              f"{np.linalg.norm(A @ x - y):.3e} {unit}")
