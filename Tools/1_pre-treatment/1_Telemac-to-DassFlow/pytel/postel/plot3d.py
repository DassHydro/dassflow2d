#!/usr/bin/python3
"""
Contains function to plot stuff
"""
from utils.exceptions import TelemacException
from postel.plot2d import set_extrema, custom_cbar
from postel.deco_cbar import deco_cbar
import numpy as np
import matplotlib.cm as cm
# Keep import enven if not used
from mpl_toolkits.mplot3d import Axes3D


def plot3d_scalar_map(fig, axe, mesh, data,
                      x_label='', y_label='', data_name='data',
                      vmin=None, vmax=None, nv=10, levels=None,
                      colorbar=True, cbar_properties=None, cbar_ticks=None,
                      cbar_ax=None, cbar_cax=None,
                      cmap_name='jet', **kwargs):
    """
    Plot a 3d representation of a triangle mesh with a data as z coordinates
    with default title (2D mesh (%d triangles, %d nodes) with "data_name" (%s))
    and xlabel and ylabel (X/Y-coordinates (m))

    @param fig (matplotlib.figure) matplotlib figure structure
    @param axe (matplotlib.axes3d) matplotlib axes on which to draw
    @param mesh (matplotlib.tri.Triangulation or tuple)
    triangular mesh or regular x,y grid
    @param data (numpy.array) Value to plot
    @param x_label (string) Label on x axe
    @param y_label (string) Label on y axe
    @param data_name (string) Name of the data to display
    @param colorbar (bool) show colorbar (default: True)
    @param cbar_ticks (list) list of values where to show color bar ticks
    @param cbar_ax (Axes) Parent axes from which space for a new colorbar
    axes will be stolen. If a list of axes is given they will all be resized
    to make room for the colorbar axes.
    @param cbar_cax (Axes) Axes into which the colorbar will be drawn.
    @param cbar_properties (dict) list additional properties of the colorbar
    @param vmin (float) Minimal value of data to plot
    @param vmax (float) Maximal value of data to plot
    @param nv (integer) Number of sample for colorbar range (default 10)
    @param levels (list) Discretisation of colorbar
    @param cmap_name (string) Name of the coloring map to use for the data
    @param kwargs (dict) rest of optional arguments given to trisurf
    """
    # Checking we indeed have a 3d Axes
    if not isinstance(axe, Axes3D):
        raise TelemacException(
                "axe must be a Axes3D instance it is %s", type(axe))

    # Getting cmap object from name
    cmap = cm.get_cmap(name=cmap_name, lut=None)

    # Setting levels with min, max and nv values if levels is not prescribed
    if vmax is not None or vmin is not None:
        assert levels is None
        vmin, vmax = set_extrema(data, vmin, vmax)
        levels = np.linspace(vmin, vmax, nv)

    # Plotting triangle surface
    img = axe.plot_trisurf(mesh, data, vmin=vmin, vmax=vmax,
                           cmap=cmap, **kwargs)

    # Adding color bar if requested
    if colorbar:
        if levels is not None and cbar_ticks is None:
            cbar_ticks = levels

        if cbar_ax is None:
            cbar_ax = axe

        if cbar_properties is not None:
            cbar = custom_cbar(
                fig, img, cbar_ax, cbar_cax, cbar_ticks, cbar_properties)
        else:
            cbar = fig.colorbar(
                img, ax=cbar_ax, cax=cbar_cax, ticks=cbar_ticks)

        if data_name is not None:
            cbar.set_label(data_name)

    # Default axis names
    axe.set_xlabel(x_label)
    axe.set_ylabel(y_label)
