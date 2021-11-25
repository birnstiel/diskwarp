import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt

# define default observer base

ex = (1, 0, 0)
ey = (0, 1, 0)
ez = (0, 0, -1)
base = np.array([ex, ey, ez])

def Rx(theta):
    nt = len(np.array(theta, ndmin=1))
    cos = np.cos(theta)
    sin = np.sin(theta)
    R = (np.eye(3)[:, :, None] * np.ones([3, 3, nt]))
    R[1, 1, :] = cos
    R[2, 2, :] = cos
    R[1, 2, :] = -sin
    R[2, 1, :] = sin
    return R
    
def Ry(theta):
    nt = len(np.array(theta, ndmin=1))
    cos = np.cos(theta)
    sin = np.sin(theta)
    R = (np.eye(3)[:, :, None] * np.ones([3, 3, nt]))
    R[0, 0, :] = cos
    R[0, 2, :] = sin
    R[2, 0, :] = -sin
    R[2, 2, :] = cos
    return R

def Rz(theta):
    nt = len(np.array(theta, ndmin=1))
    cos = np.cos(theta)
    sin = np.sin(theta)
    R = (np.eye(3)[:, :, None] * np.ones([3, 3, nt]))
    R[0, 0, :] = cos
    R[0, 1, :] = -sin
    R[1, 0, :] = sin
    R[1, 1, :] = cos
    return R

def warp(points, tilt, twist, inc=0, PA=0):
    # apply inclination
    points1 = np.einsum('ijk,klj->kli', Ry(tilt), points)
    # apply twist
    points2 = np.einsum('ijk,klj->kli', Rz(twist), points1)

    # rotate the base
    R_inc = Rx(inc)
    R_PA = Rz(PA)
    base2 = R_inc[:, :, 0].dot(R_PA[:, :, 0].dot(base))

    # project the points on observers plane
    points3 = np.einsum('ijk,kl->ijl',points2, base2)
    
    return points3

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
            _x, _y, _z =  arr[:, ir, iphi, :].T
            ax.plot(_x, _y, color=cmap(col[ir, iphi]), **kwargs)
            
    return f, ax

def get_surface(r0=0.2, r1=2, h1=0.1, hr_index=0.25, nr=20, nphi=50):


    # define cylindrical radius, azimuthal angle, and height above mid plane
    r   = np.ones([nr, nphi]) * np.linspace(r0, r1, nr)[:, None]
    phi = np.ones([nr, nphi]) * np.linspace(0, 2 * np.pi, nphi)
    z   = h1 * r**(1 + hr_index)

    # define centers as well
    rc   = 0.5 * (r[1:, 1:] + r[:-1, 1:])
    zc   = h1 * rc**(1 + hr_index)
    phic = 0.5 * (phi[1:, 1:] + phi[1:, :-1])

    # convert to cartesian (x, y = edges, xc, yc = centers)

    x   = r * np.cos(phi)
    y   = r * np.sin(phi)

    xc   = rc * np.cos(phic)
    yc   = rc * np.sin(phic)
    
    points = np.moveaxis([x, y, z], 0, 2)
    points_c = np.moveaxis([xc, yc, zc], 0, 2)
    
    return {
        'points': points,
        'points_c': points_c,
        'r': r,
        'rc': rc,
        'phi': phi,
        'phic': phic
    }
