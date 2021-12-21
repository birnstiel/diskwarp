import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt


def plot_mesh_2D(points, ax=None, **kwargs):
    """
    Takes points which is a shape (nr, nphi, 3) array
    and plots the x and y values as a grid.

    ax : plt.axes
        plot into these axes; if None, create new figure & axes
        default: None

    **kwargs are passed to the ax.plot command.

    Returns: figure, axes
    """
    if ax is None:
        f, ax = plt.subplots()
    else:
        f = ax.figure

    ax.plot(*points[:, :, :2].T, **kwargs)
    ax.plot(*np.moveaxis(points[:, :, :2], 2, 0), **kwargs)
    return f, ax


def plot_v_2D(points, v, ax=None, scale=0.05, cmap=plt.cm.RdBu_r, **kwargs):
    """
    Plot the projected velocities. The line-of-sight velocity is used as colorscale.

    points : array
        the footpoints of the vectors

    v : array
        the vectors, ankered at points

    ax : plt.axes
        plot into these axes; if None, create new figure & axes
        default: None

    **kwargs are passed to the ax.plot command.

    Returns: figure, axes
    """
    if ax is None:
        f, ax = plt.subplots()
    else:
        f = ax.figure

    arr = np.array([points, points + scale * v])
    col = Normalize()(v[:, :, -1])
    nr, nphi, _ = v.shape

    for ir in range(nr - 1):
        for iphi in range(nphi - 1):
            _x, _y, _z = arr[:, ir, iphi, :].T
            ax.plot(_x, _y, color=cmap(col[ir, iphi]), **kwargs)

    return f, ax


def get_surface(r_i, z0=0.0, psi=1.25, r_taper=80.0, q_taper=1.5, nphi=50):

    nr = len(r_i) - 1

    # define cylindrical radius, azimuthal angle, and height above mid plane
    ri = np.ones([nr + 1, nphi + 1]) * r_i[:, None]
    phii = np.ones([nr + 1, nphi + 1]) * np.linspace(0, 2 * np.pi, nphi + 1)
    zi = surface(ri, z0=z0, psi=psi, r_taper=r_taper, q_taper=q_taper)

    # define centers as well
    rc = 0.5 * (ri[1:, 1:] + ri[:-1, 1:])
    zc = surface(rc, z0=z0, psi=psi, r_taper=r_taper, q_taper=q_taper)
    phic = 0.5 * (phii[1:, 1:] + phii[1:, :-1])

    # convert to cartesian (x, y = edges, xc, yc = centers)

    xi = ri * np.cos(phii)
    yi = ri * np.sin(phii)

    xc = rc * np.cos(phic)
    yc = rc * np.sin(phic)

    points_i = np.moveaxis([xi, yi, zi], 0, 2)
    points_c = np.moveaxis([xc, yc, zc], 0, 2)

    return {
        'nr': nr,
        'nphi': nphi,
        'points_i': points_i,
        'points_c': points_c,
        'ri': ri,
        'rc': rc,
        'phii': phii,
        'phic': phic
    }

# to follow (mostly) Richs definitions


def logistic(x, a, dx, x0):
    return a / (1.0 + np.exp(-(x0 - x) / dx))


def warp(r, i_in=45.0, r0=50.0, dr=10.0):
    """return the inclination (radian) for each radius in `r`.

    Parameters
    ----------
    r : array
        radial grid
    i_in : float, optional
        maximum inner inclination in degree, by default 45.0
    r0 : float, optional
        transition radius, by default 50.0
    dr : float, optional
        transition width, by default 10.0

    Returns
    -------
    array
        inclination (radian) for each annulus in `r`
    """
    return np.radians(logistic(r, i_in, dr, r0))


def twist(r, phi=0.0, r0=50.0, dr=20.0):
    """return the twist angle (in radian) for each radius in `r`.

    Parameters
    ----------
    r : array
        radial grid
    phi : float, optional
        maximum twist of the inner disk, by default 0.0
    r0 : float, optional
        transition radius where the twist is applied, by default 50.0
    dr : float, optional
        transition width, by default 20.0

    Returns
    -------
    [type]
        [description]
    """
    return np.radians(logistic(r, phi, dr, r0))


def surface(r, z0=0.0, psi=1.25, r_taper=80.0, q_taper=1.5):
    """return a surface height `z(r)`.

    Parameters
    ----------
    r : array
        radial grid
    z0 : float, optional
        surfac height normalization at r=1, by default 0.0
    psi : float, optional
        flaring exponent, by default 1.25
    r_taper : float, optional
        radius where the surface is tapered down, by default 80.0
    q_taper : float, optional
        tapering exponent, by default 1.5

    Returns
    -------
    array
        surface height for each radius in `r`
    """
    return np.clip(z0 * r**psi * np.exp(-(r / r_taper)**q_taper), a_min=0.0, a_max=None)
