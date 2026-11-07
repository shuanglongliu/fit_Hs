"""A discretized unit sphere that turns directions into integer directory names.

Directions are snapped to a (theta, phi) grid. Both poles map to a single index (0 and the last),
so there are ``n_phi * (n_theta - 2) + 2`` directions.
"""


class Sphere:

    def __init__(self, n_theta=1801, n_phi=3600):
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.dtheta = 180 / (n_theta - 1)
        self.dphi = 360 / n_phi
        self.n_direction = n_phi * (n_theta - 2) + 2

    def index(self, direction):
        """Grid index of (theta, phi) in degrees; out-of-range angles are clamped."""
        theta, phi = direction
        theta = min(180.0, max(0.0, theta))
        phi = min(360.0, max(0.0, phi))
        i_theta = int(round(theta / self.dtheta))
        i_phi = int(round(phi / self.dphi))
        if i_phi == self.n_phi:
            i_phi = 0
        if i_theta == 0:
            return 0
        if i_theta == self.n_theta - 1:
            return self.n_direction - 1
        return self.n_phi * (i_theta - 1) + i_phi + 1

    def direction(self, index):
        """(theta, phi) in degrees of a grid index; out-of-range indices are clamped."""
        index = min(self.n_direction - 1, max(0, index))
        if index == 0:
            return (0.0, 0.0)
        if index == self.n_direction - 1:
            return (180.0, 0.0)
        i_theta = (index - 1) // self.n_phi + 1
        i_phi = (index - 1) % self.n_phi
        return (i_theta * self.dtheta, i_phi * self.dphi)
