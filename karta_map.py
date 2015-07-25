""" Karta-map is a module that depends on the karta package and adds some
useful mapping functions. """

import scipy.optimize
import karta
import numpy as np

from typing import Iterable
from matplotlib import patches
from matplotlib.pyplot import gca, Axes

def get_axes_extents(ax, ax_crs: karta.crs.CRS, crs=karta.crs.SphericalEarth):
    """ Get the extents of an Axes in geographical (or other) coordinates. """
    xl, xr = ax.get_xlim()
    yb, yt = ax.get_ylim()
    
    transform = lambda x, y: crs.project(*ax_crs.project(x, y, inverse=True))
    ll = transform(xl, yb)
    lr = transform(xr, yb)
    ur = transform(xr, yt)
    ul = transform(xl, yt)
    return karta.Polygon([ll, lr, ur, ul], crs=crs)

def _segment(ax,
             x0: float, 
             x1: float,
             y0: float,
             y1: float,
             n: int,
             crs: karta.crs.CRS, **kw):
    """ Convert a line between two points in geographical space into a
    (segmented) curve between those points in a projected space. """
    x = np.linspace(x0, x1, n)
    y = np.linspace(y0, y1, n)
    xp, yp = crs.project(x, y)
    return ax.plot(xp, yp, **kw)

def curve(ax, x, y, n=20, crs=karta.crs.SphericalEarth, **kw):
    segs = []
    for i in range(len(x)-1):
        segs.append(_segment(ax, x[i], x[i+1], y[i], y[i+1], n, crs, **kw))
    return segs

def add_graticule(ax, xs: Iterable, ys: Iterable,
                  ax_crs=karta.crs.Cartesian,
                  graticule_crs=karta.crs.SphericalEarth,
                  lineargs=None):

    if lineargs is None:
        lineargs = dict(color="k", linewidth=0.5)

    bbox = get_axes_extents(ax, ax_crs, graticule_crs)
    x, y = bbox.get_coordinate_lists(crs=graticule_crs)
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    for i in range(len(xs)):
        curve(ax, (xs[i], xs[i]), (ymin, ymax), crs=ax_crs, **lineargs)
    for i in range(len(ys)):
        curve(ax, (xmin, xmax), (ys[i], ys[i]), crs=ax_crs, **lineargs)
    return

def find_intersection(xa, ya, xb, yb, crsa, crsb):
    """ Return the point in *crs_a* on the line xa, ya that intersects xb, yb on *crs_b*. """


def isbetween(x: float, a: float, b: float) -> bool:
    return (a < x < b) or (b < x < a)

def froot(f, a, b):
    return scipy.optimize.brentq(f, a, b)

def label_ticks(ax, xs: Iterable, ys: Iterable,
                ax_crs=karta.crs.Cartesian,
                graticule_crs=karta.crs.SphericalEarth,
                textargs=None, tickargs=None,
                x_suffix="\u00b0E", y_suffix="\u00b0N"):

    if textargs is None:
        textargs = dict()

    if tickargs is None:
        tickargs = dict(marker="+", mew=2, ms=14, mfc="k", mec="k", ls="none")

    # Find tick locations
    bbox = get_axes_extents(ax, ax_crs, graticule_crs)  # bottom, right, top, left

    ticks = dict(xticks=[], yticks=[])

    xmin, xmax = sorted(ax.get_xlim())
    ymin, ymax = sorted(ax.get_ylim())

    tickproj = graticule_crs.project
    axproj = ax_crs.project
    #ax_inv = lambda x, y: ax_crs.project(x, y, inverse=True)

    # bottom spine
    for x in xs:
        if isbetween(x, bbox[0][0], bbox[1][0]):
            ticks["xticks"].append((froot(lambda xt: tickproj(*axproj(xt, ymin, inverse=True))[0]-x, xmin, xmax),
                                    ymin,
                                    "{0}{1}".format(x, x_suffix)))

    for y in ys:
        if isbetween(y, bbox[0][1], bbox[1][1]):
            ticks["yticks"].append((froot(lambda xt: tickproj(*axproj(xt, ymin, inverse=True))[1]-y, xmin, xmax),
                                    ymin,
                                    "{0}{1}".format(y, y_suffix)))

    # top spine
    for x in xs:
        if isbetween(x, bbox[2][0], bbox[3][0]):
            ticks["xticks"].append((froot(lambda xt: tickproj(*axproj(xt, ymax, inverse=True))[0]-x, xmin, xmax),
                                    ymax,
                                    "{0}{1}".format(x, x_suffix)))

    for y in ys:
        if isbetween(y, bbox[2][1], bbox[3][1]):
            ticks["yticks"].append((froot(lambda xt: tickproj(*axproj(xt, ymax, inverse=True))[1]-y, xmin, xmax),
                                    ymax,
                                    "{0}{1}".format(y, y_suffix)))

    # left spine
    for x in xs:
        if isbetween(x, bbox[0][0], bbox[3][0]):
            ticks["xticks"].append((xmin,
                                    froot(lambda yt: tickproj(*axproj(xmin, yt, inverse=True))[0]-x, ymin, ymax),
                                    "{0}{1}".format(x, x_suffix)))


    for y in ys:
        if isbetween(y, bbox[0][1], bbox[3][1]):
            ticks["yticks"].append((xmin,
                                    froot(lambda yt: tickproj(*axproj(xmin, yt, inverse=True))[1]-y, ymin, ymax),
                                    "{0}{1}".format(y, y_suffix)))


    # right spine
    for x in xs:
        if isbetween(x, bbox[1][0], bbox[2][0]):
            ticks["xticks"].append((xmax,
                                    froot(lambda yt: tickproj(*axproj(xmax, yt, inverse=True))[0]-x, ymin, ymax),
                                    "{0}{1}".format(x, x_suffix)))

    for y in ys:
        if isbetween(y, bbox[1][1], bbox[2][1]):
            ticks["yticks"].append((xmax,
                                    froot(lambda yt: tickproj(*axproj(xmax, yt, inverse=True))[1]-y, ymin, ymax),
                                    "{0}{1}".format(y, y_suffix)))

    # Update map

    for pt in ticks["xticks"]:
        ax.plot(pt[0], pt[1], **tickargs)
        ax.text(pt[0], pt[1], pt[2], **textargs)

    for pt in ticks["yticks"]:
        ax.plot(pt[0], pt[1], **tickargs)
        ax.text(pt[0], pt[1], pt[2], **textargs)

    # # Apply to the Axes
    ax.set_xticks([])
    ax.set_yticks([])

    # ax.set_xticklabels(xs)
    # ax.set_yticklabels(ys)
    return

def plot(geoms: Iterable, *args, ax=None, crs=None, **kwargs):
    """ Metafunction that dispatches to the correct plotting routine. """

    ax = gca()

    if isinstance(geoms, list):
        if geoms[0]._geotype == "Point":
            ret = plot_points(geoms, *args, ax=ax, crs=crs, **kwargs)
        elif geoms[0]._geotype == "Multipoint":
            ret = plot_multipoints(geoms, *args, ax=ax, crs=crs, **kwargs)
        elif geoms[0]._geotype == "Line":
            ret = plot_lines(geoms, *args, ax=ax, crs=crs, **kwargs)
        elif geoms[0]._geotype == "Polygon":
            ret = plot_polygons(geoms, ax=ax, crs=crs, **kwargs)
        else:
            raise TypeError("Invalid geotype")
    else:
        if geoms._geotype == "Point":
            ret = plot_point(geoms, *args, ax=ax, crs=crs, **kwargs)
        elif geoms._geotype == "Multipoint":
            ret = plot_multipoint(geoms, *args, ax=ax, crs=crs, **kwargs)
        elif geoms._geotype == "Line":
            ret = plot_line(geoms, *args, ax=ax, crs=crs, **kwargs)
        elif geoms._geotype == "Polygon":
            ret = plot_polygon(geoms, ax=ax, crs=crs, **kwargs)
        else:
            raise TypeError("Invalid geotype")

    return ret

def plot_line(geom, ax, *args, crs=None, **kwargs):
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

def plot_lines(geoms, ax, *args, crs=None, **kwargs):
    return [plot_line(geom, ax, *args, crs=crs, **kwargs) for geom in geoms]

def plot_polygon(geom, ax, crs=None, **kwargs):
    p = patches.Polygon(geom.get_vertices(crs=crs), **kwargs)
    ax.add_patch(p)
    return p

def plot_polygons(geoms: Iterable, ax, crs=None, **kwargs):
    """ Plot Polygon geometries, projected to the coordinate system `crs` """
    return [plot_polygon(geom, ax, crs=crs, **kwargs) for geom in geoms]

