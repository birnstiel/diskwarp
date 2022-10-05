import subprocess

import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from tqdm.auto import tqdm


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


def logistic(a, r, r0, dr):
    """Logistic function with steeper transition

    Parameters
    ----------
    a : float
        maximum warp angle in DEGREE
    r : array
        radial grid
    r0 : float
        transition radius
    dr : float
        transition width. after this width, we are within
        1/(1+exp(10)) < 5e-5 of the final radius

    Returns
    -------
    array
        the inclination array in RADIAN for every radius in `r`
    """
    return np.radians(a / (1 + np.exp((r - r0) / (0.1 * dr))))


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


def warp_transformation(x, y, z, theta, phi):
    """applies a warp to the cartesian coordinates `x, y, z` by warping and twisting by theta and phi (in RADIAN)

    Parameters
    ----------
    x : array
        x-coordinates
    y : array
        y-coorindates
    z : array
        z-coorindates
    theta : float | array
        the inclination angle of the warp in (RADIAN)
    phi : float | array
        the twist angle of the warp (in RADIAN)

    Returns
    -------
    x,y,z
        tuple with the updated cartesian coordinates
    """
    xprime = x * np.cos(phi) - y * np.sin(phi) * np.cos(theta) + z * np.sin(phi) * np.sin(theta)
    yprime = x * np.sin(phi) + y * np.cos(phi) * np.cos(theta) - z * np.sin(theta) * np.cos(phi)
    zprime = y * np.sin(theta) + z * np.cos(theta)
    return xprime, yprime, zprime


def unwarp_transformation(x, y, z, theta, phi):
    """undoes a warp to the cartesian coordinates `x, y, z` by un-warping and un-twisting
    a warp/twist with angles theta and phi (in RADIAN)

     Parameters
     ----------
     x : array
         x-coordinates
     y : array
         y-coorindates
     z : array
         z-coorindates
     theta : float | array
         the inclination angle of the warp in (RADIAN)
     phi : float | array
         the twist angle of the warp (in RADIAN)

     Returns
     -------
     x,y,z
         tuple with the updated cartesian coordinates
     """
    xprime = x * np.cos(phi) + y * np.sin(phi)
    yprime = -x * np.sin(phi) * np.cos(theta) + y * np.cos(phi) * np.cos(theta) + z * np.sin(theta)
    zprime = x * np.sin(phi) * np.sin(theta) - y * np.sin(theta) * np.cos(phi) + z * np.cos(theta)
    return xprime, yprime, zprime


def vel_sph_to_car(theta, phi, vr, vtheta, vphi):
    """given theta, phi, convert spherical velocities to cartesian,

    Parameters
    ----------
    theta : array
        theta grid
    phi : array
        phi grid
    vr, vtheta, vphi : arrays
        radial, theta, and phi components of the velocities

    Returns
    -------
    vx, vy, vz
        cartesian x-, y-, and z-components of the velocities
    """
    vx = vr * np.sin(theta) * np.cos(phi) + vtheta * np.cos(phi) * np.cos(theta) - vphi * np.sin(phi)
    vy = vr * np.sin(phi) * np.sin(theta) + vtheta * np.sin(phi) * np.cos(theta) + vphi * np.cos(phi)
    vz = vr * np.cos(theta) - vtheta * np.sin(theta)

    return vx, vy, vz


def vel_car_to_sph(theta, phi, vx, vy, vz):
    """given theta, phi, convert cartesian velocities to spherical

    Parameters
    ----------
    theta : array
        theta grid
    phi : array
        phi grid
    vx, vy, vz : arrays
        cartesian x-, y-, and z-components of the velocities

    Returns
    -------
    vr, vtheta, vphi
        spherical r-, theta-, and phi-components of the velocities
    """
    vr = vx * np.sin(theta) * np.cos(phi) + vy * np.sin(theta) * np.sin(phi) + vz * np.cos(theta)
    vt = vx * np.cos(theta) * np.cos(phi) + vy * np.cos(theta) * np.sin(phi) - vz * np.sin(theta)
    vp = -vx * np.sin(phi) + vy * np.cos(phi)

    return vr, vt, vp


def call_radmc(cmd, verbose=False, total=None):
    """
    Run radmc3d command and show progress bar instead.

    cmd : str
        the command to run, e.g. 'radmc3d mctherm'

    verbose : bool
        if True, then all output except the photon packges are shown
        if False, just the progress is shown.

    total : None | int
        total number of photon packages, if known
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    output = []

    if 'nphot' in cmd:
        total = int(cmd.split('nphot')[1].split()[1])

    with tqdm(total=total, unit='photons') as pbar:
        for line in p.stdout:
            if 'Photon nr' in line:
                pbar.update(1000)
            elif verbose:
                print(line, end='')
            output += [line]
    rc = p.wait()
    return rc, ''.join(output)


def grid_refine_inner_edge(x_orig, nlev, nspan):
    "subdivide the `nspan` cells in grid `x_orig` by `nlev` levels"
    x = x_orig.copy()
    rev = x[0] > x[1]
    for ilev in range(nlev):
        x_new = 0.5 * (x[1:nspan + 1] + x[:nspan])
        x_ref = np.hstack((x, x_new))
        x_ref.sort()
        x = x_ref
        if rev:
            x = x[::-1]
    return x
