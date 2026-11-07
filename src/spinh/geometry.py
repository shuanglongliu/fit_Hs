"""Spherical coordinates, reference frames and rotations.

Angles are in degrees unless ``rad=True``:
    theta : polar angle from +z, in [0, 180]
    phi   : azimuth in the xy plane from +x toward +y, in [0, 360)

A reference frame is a 3x3 array ``emat`` whose rows are the frame axes (ex, ey, ez) in global
coordinates, so the components of a global vector v in that frame are ``emat @ v``.
"""

import math

import numpy as np
from scipy.spatial.transform import Rotation


def sph2cart(sph, rad=False):
    """[r, theta, phi] -> np.array([x, y, z])."""
    r, theta, phi = sph
    if not rad:
        theta, phi = np.deg2rad(theta), np.deg2rad(phi)
    return np.array([r * math.sin(theta) * math.cos(phi),
                     r * math.sin(theta) * math.sin(phi),
                     r * math.cos(theta)])


def cart2sph(xyz, rad=False):
    """[x, y, z] -> [r, theta, phi], with phi in (-180, 180]."""
    xyz = np.array([xyz], dtype=float)
    xy = xyz[:, 0]**2 + xyz[:, 1]**2
    r = np.sqrt(xy + xyz[:, 2]**2)[0]
    theta = np.arctan2(np.sqrt(xy), xyz[:, 2])[0]
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])[0]
    if rad:
        return [r, theta, phi]
    return [r, theta / math.pi * 180, phi / math.pi * 180]


def opposite_direction(direction):
    """(theta, phi) of the opposite direction, phi in [0, 360) and 0 at the poles."""
    theta, phi = direction
    theta = 180 - theta
    phi = (phi + 180) % 360
    if theta == 0 or theta == 180:
        phi = 0
    return (theta, phi)


def angle_between(a, b, rad=False):
    """Angle between two cartesian vectors."""
    cosab = np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b)
    angle = math.acos(min(1, max(-1, cosab)))
    return angle if rad else 180 * angle / math.pi


def angle_between_sph(a, b, rad=False):
    """Angle between two vectors given as [r, theta, phi]."""
    return angle_between(sph2cart(a, rad), sph2cart(b, rad), rad)


def change_frame_sph(vi_sph, emati, ematf, rad=False):
    """Re-express [r, theta, phi] given in frame ``emati`` in frame ``ematf``.

    phi of the result is in [0, 360) (or [0, 2 pi)) and 0 at the poles.
    """
    v = np.matmul(np.transpose(emati), sph2cart(vi_sph, rad=rad))
    v = np.matmul(np.linalg.inv(np.transpose(ematf)), v)
    vf_sph = cart2sph(v, rad=rad)

    pi, twopi = (np.pi, 2 * np.pi) if rad else (180, 360)
    if abs(vf_sph[1]) < 1.e-9 or abs(vf_sph[1] - pi) < 1.e-9:
        vf_sph[2] = 0
    if vf_sph[2] < 0:
        vf_sph[2] = vf_sph[2] + twopi
    elif vf_sph[2] > twopi:
        vf_sph[2] = vf_sph[2] - twopi
    return vf_sph


def change_frame_angles(angles, emati, ematf):
    """Re-express arrays (thetas, phis), in radians, from frame ``emati`` in frame ``ematf``."""
    thetas, phis = angles
    new = [change_frame_sph([1, t, p], emati, ematf, rad=True)[1:] for t, p in zip(thetas, phis)]
    new = np.array(new).reshape(-1, 2)
    return new[:, 0], new[:, 1]


def get_emat_local(z_direction=(0, 0)):
    """Right-handed frame whose z axis points along ``z_direction = (theta, phi)``.

    The global x axis lies in the local xz plane. ``z_direction == (90, 0)`` (a tuple) is the
    degenerate case handled explicitly: local (x, y, z) = global (y, z, x).
    """
    if z_direction == (90, 0):
        return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    ez = sph2cart([1, z_direction[0], z_direction[1]])
    ey = np.cross(ez, [1, 0, 0])
    ey = ey / np.linalg.norm(ey)
    ex = np.cross(ey, ez)
    return np.array([ex, ey, ez])


def euler_zyz(alpha, beta, gamma, degrees=True):
    """Rotation matrix Rz(gamma) Ry(beta) Rz(alpha): extrinsic zyz Euler angles.

    Applied to the z axis it gives the direction (theta, phi) = (beta, gamma); alpha is the turn
    about z made first.
    """
    return Rotation.from_euler("zyz", [alpha, beta, gamma], degrees=degrees).as_matrix()


def rotate_frame(emat, alpha, beta, gamma):
    """Rotate the frame ``emat`` by the extrinsic zyz Euler angles (radians) about its own axes.

    The new z axis points along (theta, phi) = (beta, gamma) in the old frame, and alpha turns the
    new x and y axes about it.
    """
    return euler_zyz(alpha, beta, gamma, degrees=False).T @ np.asarray(emat, dtype=float)
