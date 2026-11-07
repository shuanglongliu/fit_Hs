"""Ready-made directions and configurations for the noncollinear job sets.

Directions are (theta, phi) in degrees. ``meridian``, ``parallel`` and ``great_circles`` return
lists of directions, for ``jobs.add_uniform`` (every site along each direction); ``updown``
returns configurations, for ``jobs.add``.
"""

import numpy as np

X, Y, Z = (90.0, 0.0), (90.0, 90.0), (0.0, 0.0)
UP, DOWN = (0.0, 0.0), (180.0, 0.0)


def meridian(phi, theta_min=0, theta_max=180, n=7):
    """n directions from theta_min to theta_max (both included) at fixed phi."""
    return [(theta, phi) for theta in np.linspace(theta_min, theta_max, n)]


def parallel(theta, phi_min=0, phi_max=330, n=12):
    """n directions from phi_min to phi_max (both included) at fixed theta."""
    return [(theta, phi) for phi in np.linspace(phi_min, phi_max, n)]


def great_circles(n_theta=7, n_phi=12):
    """The meridians phi = 0, 180, 90, 270 (n_theta points each, poles included) and the equator
    (n_phi points); repeated directions are dropped by ``jobs.add_uniform``."""
    directions = []
    for phi in (0, 180, 90, 270):
        directions += meridian(phi, 0, 180, n_theta)
    return directions + parallel(90, 0, 360 - 360 / n_phi, n_phi)


def updown(configs):
    """Collinear configurations from strings of 'u' (theta = 0) and 'd' (theta = 180)."""
    for config in configs:
        if set(config) - set("ud"):
            raise ValueError(f"Configuration {config!r}: use only 'u' and 'd'.")
    return [[UP if c == "u" else DOWN for c in config] for config in configs]
