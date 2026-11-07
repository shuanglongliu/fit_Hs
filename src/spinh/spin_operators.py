"""Matrices of spin operators in the |S, m> basis, m = S, S-1, ..., -S."""

from types import SimpleNamespace

import numpy as np


def spin_operators(S):
    """Sx, Sy, Sz (hbar = 1), the identity ID, dim = 2S + 1 and ss = S(S + 1)."""
    dim = int(round(2 * S + 1))
    m = S - np.arange(dim)
    Sz = np.diag(m).astype(complex)
    # <m + 1| S+ |m> = sqrt(S(S + 1) - m(m + 1))
    Splus = np.diag(np.sqrt(S * (S + 1) - m[1:] * (m[1:] + 1)), k=1).astype(complex)
    Sminus = Splus.conj().T
    Sx = (Splus + Sminus) / 2
    Sy = (Splus - Sminus) / 2j
    return SimpleNamespace(S=S, dim=dim, ss=S * (S + 1), ID=np.eye(dim, dtype=complex),
                           Sx=Sx, Sy=Sy, Sz=Sz, Splus=Splus, Sminus=Sminus)


def top_state(op):
    """Eigenvector of the Hermitian matrix ``op`` with the largest eigenvalue."""
    w, v = np.linalg.eigh(op)
    return v[:, np.argmax(w)]


def expectation(state, op):
    return np.real(np.vdot(state, op @ state))
