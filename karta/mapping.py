"""
karta.mapping is an extension module for the karta package that provides a
(thin) wrapper around matplotlib for easily plotting geogrphaical data.
"""

import scipy.optimize
import numpy as np

from typing import Iterable
from matplotlib import patches
from matplotlib.pyplot import gca, Axes

from .vector import Point, Multipoint, Line, Polygon, Geometry
from .raster import RegularGrid
from .crs import CRS, Cartesian, SphericalEarth

def default_current_axes(wrappedfunc):
    """ Decorator to set current Axes as default ax in plotting functions """
    def replacementfunc(*args, **kwargs):
        kwargs.setdefault("ax", gca())
        return wrappedfunc(*args, **kwargs)
    return replacementfunc

def recurse_iterables(wrappedfunc):
    """ Decorator to generate functions that apply themselve recursively to
    iterable non-Geometry inputs. """
    def replacementfunc(main_arg, *args, **kwargs):
        if not isinstance(main_arg, Geometry) and hasattr(main_arg, "__iter__"):
            return [replacementfunc(a, *args, **kwargs) for a in main_arg]
        else:
            return wrappedfunc(main_arg, *args, **kwargs)
    return replacementfunc

def get_axes_extent(ax, ax_crs: CRS, crs=SphericalEarth):
    """ Get the extent of an Axes in geographical (or other) coordinates. """
    xl, xr = ax.get_xlim()
    yb, yt = ax.get_ylim()

    ll = ax_crs.transform(crs, xl, yb)
    lr = ax_crs.transform(crs, xr, yb)
    ur = ax_crs.transform(crs, xr, yt)
    ul = ax_crs.transform(crs, xl, yt)
    return Polygon([ll, lr, ur, ul], crs=crs)

def geodesic(pt0: Point, pt1: Point, n=20):
    """ Return a Line representing the geodesic path between two points, approximated by `n` segments. """
    i = 0
    cur_pt = pt0
    points = [cur_pt]
    while i != n:
        remaining_dist = cur_pt.distance(pt1)
        step_dist = remaining_dist / (n-i)
        az = cur_pt.azimuth(pt1)
        cur_pt = cur_pt.walk(step_dist, az)
        points.append(cur_pt)
        i += 1
    return Line(points)

def add_graticule(ax, xs: Iterable, ys: Iterable,
                  map_crs=Cartesian,
                  graticule_crs=SphericalEarth,
                  lineargs=None):
    """ Add a map graticule, with intervals in `graticule_crs` projected onto a
    map projected with `map_crs` """
    if lineargs is None:
        lineargs = dict(color="k", linewidth=0.5)

    bbox = get_axes_extent(ax, map_crs, graticule_crs)
    x, y = bbox.get_coordinate_lists(crs=graticule_crs)
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    for i in range(len(xs)):
        _line = Line([(xs[i], ymin), (xs[i], ymax)], crs=map_crs)
        plot(_line, ax, **lineargs)
    for i in range(len(ys)):
        _line = Line([(xmin, ys[i]), (xmax, ys[i])], crs=map_crs)
        plot(_line, ax, **lineargs)
    return

def add_graticule_contour(ax, xs, ys, map_crs, graticule_crs, nx=100, ny=100):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmap = np.linspace(xmin, xmax, nx)
    ymap = np.linspace(ymin, ymax, ny)
    Xm, Ym = np.meshgrid(xmap, ymap)
    Xg, Yg = map_crs.transform(graticule_crs, Ym, Ym)
    ax.contour(Xm, Ym, abs(Xg), levels=xs, colors="k", linestyles="-")
    ax.contour(Xm, Ym, Yg, levels=ys, colors="k", linestyles="-")

    Xg_pm = Xg
    Xg_pm[(abs(Xg)>10) & (abs(Xg)<170)] = np.nan
    ax.contour(Xm, Ym, Xg_pm, levels=[0.0], colors="k", linestyles="-")

    return

def find_intersection(xa, ya, xb, yb, crsa, crsb):
    """ Return the point in *crs_a* on the line xa, ya that intersects xb, yb on *crs_b*. """


def isbetween(x: float, a: float, b: float) -> bool:
    return (a < x < b) or (b < x < a)

def froot(f, a, b):
    return scipy.optimize.brentq(f, a, b)

def label_ticks(ax, xs: Iterable, ys: Iterable,
                map_crs=Cartesian,
                graticule_crs=SphericalEarth,
                textargs=None, tickargs=None,
                x_suffix="\u00b0E", y_suffix="\u00b0N"):

    if textargs is None:
        textargs = dict()

    if tickargs is None:
        tickargs = dict(marker="+", mew=2, ms=14, mfc="k", mec="k", ls="none")

    # Find tick locations
    bbox = get_axes_extent(ax, map_crs, graticule_crs)  # bottom, right, top, left

    ticks = dict(xticks=[], yticks=[])

    xmin, xmax = sorted(ax.get_xlim())
    ymin, ymax = sorted(ax.get_ylim())

    # bottom spine
    for x in xs:
        if isbetween(x, bbox[0][0], bbox[1][0]):
            ticks["xticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymin)[0]-x,
                                          xmin, xmax),
                                    ymin,
                                    "{0}{1}".format(x, x_suffix)))

    for y in ys:
        if isbetween(y, bbox[0][1], bbox[1][1]):
            ticks["yticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymin)[1]-y,
                                          xmin, xmax),
                                    ymin,
                                    "{0}{1}".format(y, y_suffix)))

    # top spine
    for x in xs:
        if isbetween(x, bbox[2][0], bbox[3][0]):
            ticks["xticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymax)[0]-y,
                                          xmin, xmax),
                                    ymax,
                                    "{0}{1}".format(x, x_suffix)))

    for y in ys:
        if isbetween(y, bbox[2][1], bbox[3][1]):
            ticks["yticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymax)[1]-y,
                                          xmin, xmax),
                                    ymax,
                                    "{0}{1}".format(y, y_suffix)))

    # left spine
    for x in xs:
        if isbetween(x, bbox[0][0], bbox[3][0]):
            ticks["xticks"].append((xmin,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmin, yt)[0]-x,
                                          ymin, ymax),
                                    "{0}{1}".format(x, x_suffix)))


    for y in ys:
        if isbetween(y, bbox[0][1], bbox[3][1]):
            ticks["yticks"].append((xmin,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmin, yt)[1]-y,
                                          ymin, ymax),
                                    "{0}{1}".format(y, y_suffix)))


    # right spine
    for x in xs:
        if isbetween(x, bbox[1][0], bbox[2][0]):
            ticks["xticks"].append((xmax,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmax, yt)[0]-x,
                                          ymin, ymax),
                                    "{0}{1}".format(x, x_suffix)))

    for y in ys:
        if isbetween(y, bbox[1][1], bbox[2][1]):
            ticks["yticks"].append((xmax,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmax, yt)[1]-y,
                                          ymin, ymax),
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
    return

def _get_plotting_func(geom):
    if isinstance(geom, list):
        return _get_plotting_func(geom[0])
    if not hasattr(geom, "_geotype"):
        raise TypeError("Invalid input type: {0}".format(type(geom)))
    if geom._geotype == "Point":
        return plot_point
    if geom._geotype == "Multipoint":
        return plot_multipoint
    if geom._geotype == "Line":
        return plot_line
    if geom._geotype == "Polygon":
        return plot_polygon
    raise TypeError("Invalid geotype: {0}".format(geom._geotype))

def plot(geom, *args, **kwargs):
    """ Metafunction that dispatches to the correct plotting routine. """
    func = _get_plotting_func(geom)
    return func(geom, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_point(geom, *args, ax=None, crs=None, **kwargs):
    """ Plot a Point geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("marker", ".")
    x, y = geom.get_vertex(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_multipoint(geom, *args, ax=None, crs=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("linestyle", "none")
    kwargs.setdefault("marker", ".")
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_line(geom, *args, ax=None, crs=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_polygon(geom, *args, ax=None, crs=None, **kwargs):
    """ Plot a Polygon geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("facecolor", "none")
    kwargs.setdefault("edgecolor", "black")
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.fill(x, y, *args, **kwargs)

@default_current_axes
def plot_grid(grid: RegularGrid, ax=None, crs=None, **kwargs):
    kwargs.setdefault("origin", "bottom")
    kwargs.setdefault("extent", grid.get_extent(crs=crs))
    return ax.imshow(grid.values, **kwargs)

def _position_over(artist):
    xy = artist.get_xy()
    x = xy[:,0]
    y = xy[:,1]
    return 0.5*(np.min(x) + np.max(x)), 0.5*(np.min(y) + np.max(y))

def _position_below(artist):
    xy = artist.get_xy()
    x = xy[:,0]
    y = xy[:,1]
    return 0.5*(np.min(x) + np.max(x)), np.min(y)

def _position_above(artist):
    xy = artist.get_xy()
    x = xy[:,0]
    y = xy[:,1]
    return 0.5*(np.min(x) + np.max(x)), np.max(y)

@default_current_axes
def annotate(artist, label, where="over", ax=None, **kwargs):
    """ Add a Text object near *artist*. """
    if where == "over":
        x, y = _position_over(artist)
        kwargs.setdefault("va", "center")
    elif where == "below":
        x, y = _position_below(artist)
        kwargs.setdefault("va", "top")
    elif where == "above":
        x, y = _position_above(artist)
        kwargs.setdefault("va", "bottom")
    else:
        raise ValueError("invalid value for 'where'")
    return ax.text(x, y, label, **kwargs)

