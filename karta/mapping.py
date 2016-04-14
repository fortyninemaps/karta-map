"""
karta.mapping is an extension module for the karta package that provides a
(thin) wrapper around matplotlib for easily plotting geogrphaical data.
"""

import scipy.optimize
import numpy as np

from typing import Union, Iterable, Tuple, Callable
from matplotlib.pyplot import gca, Axes, Artist

from .vector import Point, Multipoint, Line, Polygon, Geometry
from .raster import RegularGrid
from .crs import CRS, Cartesian, SphericalEarth

def default_current_axes(wrappedfunc: Callable):
    """ Decorator to set current Axes as default ax in plotting functions """
    def replacementfunc(*args, **kwargs):
        if "ax" not in kwargs:
            # don't use dict.setdefault because don't want gca() to be evaluated
            kwargs["ax"] = gca()
        return wrappedfunc(*args, **kwargs)
    return replacementfunc

def recurse_iterables(wrappedfunc: Callable):
    """ Decorator to generate functions that apply themselve recursively to
    iterable non-Geometry inputs. """
    def replacementfunc(main_arg, *args, **kwargs):
        if not isinstance(main_arg, Geometry) and hasattr(main_arg, "__iter__"):
            return [replacementfunc(a, *args, **kwargs) for a in main_arg]
        else:
            return wrappedfunc(main_arg, *args, **kwargs)
    return replacementfunc

def get_axes_extent(ax: Axes, ax_crs: CRS, crs: CRS=SphericalEarth):
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

@default_current_axes
def add_graticule(xs: Iterable[float], ys: Iterable[float],
        ax: Axes=None, map_crs: CRS=Cartesian, graticule_crs: CRS=SphericalEarth,
        lineargs=None):
    """ Adds a map graticule.

    Parameters
    ----------
    xs : Iterable[float],
    ys : Iterable[float]
        Easting and northing componenets of graticule, in `graticule_crs`
    ax : Axes, optional
        Axes to draw to (default current Axes)
    map_crs : karta.crs.CRS, optional
        CRS defining the display projection (default Cartesian)
    graticule_crs : karta.crs.CRS, optional
        CRS defining the graticule projection (default SphericalEarth)
    lineargs : dict, optional
        Arguments passed to karta.mapping.plot while drawing graticule lines
    """
    if lineargs is None:
        lineargs = dict(color="k", linewidth=0.5)

    bbox = get_axes_extent(ax, map_crs, graticule_crs)
    x, y = bbox.get_coordinate_lists(crs=graticule_crs)
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    artists = []
    for i in range(len(xs)):
        _line = Line([(xs[i], ymin), (xs[i], ymax)], crs=graticule_crs)
        artists.append(plot(_line, ax=ax, crs=map_crs, **lineargs))
    for i in range(len(ys)):
        _line = Line([(xmin, ys[i]), (xmax, ys[i])], crs=graticule_crs)
        artists.append(plot(_line, ax=ax, crs=map_crs, **lineargs))
    return artists

def isbetween(x: float, a: float, b: float) -> bool:
    return (a < x < b) or (b < x < a)

def froot(f: float, a: float, b: float) -> float:
    return scipy.optimize.brentq(f, a, b)

@default_current_axes
def label_ticks(xs: Iterable[float], ys: Iterable[float],
        ax: Axes=None, map_crs: CRS=Cartesian, graticule_crs: CRS=SphericalEarth,
        textargs=None, tickargs=None,
        x_suffix: str="E", y_suffix: str="N"):
    """ Label graticule lines, returning a list if Text objects.

    Parameters
    ----------
    xs : Iterable[float],
    ys : Iterable[float]
        Easting and northing componenets of labels, in `graticule_crs`
    ax : Axes, optional
        Axes to draw to (default current Axes)
    map_crs : karta.crs.CRS, optional
        CRS giving the display projection (default Cartesian)
    graticule_crs : karta.crs.CRS, optional
        CRS giving the graticule/label projection (default SphericalEarth)
    textargs : dict, optional
        Keyword arguments to pass to plt.text
    tickargs : dict, optional
        Keyword arguments to pass to plt.plot
    x_suffix : str, optional
        Suffix for eastings labels (default 'E')
    y_suffix : str, optional
        Suffix to pass to northings labels (default 'N')
    """

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
                                          map_crs.transform(graticule_crs, xt, ymax)[0]-x,
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
    txts = []
    for pt in ticks["xticks"]:
        ax.plot(pt[0], pt[1], **tickargs)
        txts.append(ax.text(pt[0], pt[1], pt[2], **textargs))

    for pt in ticks["yticks"]:
        ax.plot(pt[0], pt[1], **tickargs)
        txts.append(ax.text(pt[0], pt[1], pt[2], **textargs))

    ax.set_xticks([])
    ax.set_yticks([])
    return txts

def _get_plotting_func(geom: Union[Geometry, Iterable[Geometry]]) -> Callable:
    if isinstance(geom, RegularGrid):
        return plot_grid
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

def plot(geom: Union[Union[Geometry, RegularGrid], Iterable[Union[Geometry, RegularGrid]]], *args, **kwargs):
    """ Metafunction that dispatches to the correct plotting routine. """
    func = _get_plotting_func(geom)
    return func(geom, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_point(geom: Union[Point, Iterable[Point]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Point geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("marker", ".")
    x, y = geom.get_vertex(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_multipoint(geom: Union[Multipoint, Iterable[Multipoint]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("linestyle", "none")
    kwargs.setdefault("marker", ".")
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_line(geom: Union[Line, Iterable[Line]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
@recurse_iterables
def plot_polygon(geom: Union[Polygon, Iterable[Polygon]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Polygon geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("facecolor", "none")
    kwargs.setdefault("edgecolor", "black")
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.fill(x, y, *args, **kwargs)

@default_current_axes
def plot_grid(grid: RegularGrid, ax: Axes=None, crs: CRS=None, **kwargs):
    kwargs.setdefault("origin", "bottom")
    kwargs.setdefault("extent", grid.get_extent(crs=crs))

    # compute the pixels that can actually be displayed
    _, _, width, height = ax.bbox.bounds
    ny, nx = grid.size
    r = (max(int(ny//height), 1), max(int(nx//width), 1))
    return ax.imshow(grid[::r[0],::r[1]], **kwargs)

def _position_over(artist: Artist) -> Tuple[float, float]:
    xy = artist.get_xy()
    x = xy[:,0]
    y = xy[:,1]
    return 0.5*(np.min(x) + np.max(x)), 0.5*(np.min(y) + np.max(y))

def _position_below(artist: Artist) -> Tuple[float, float]:
    xy = artist.get_xy()
    x = xy[:,0]
    y = xy[:,1]
    return 0.5*(np.min(x) + np.max(x)), np.min(y)

def _position_above(artist: Artist) -> Tuple[float, float]:
    xy = artist.get_xy()
    x = xy[:,0]
    y = xy[:,1]
    return 0.5*(np.min(x) + np.max(x)), np.max(y)

@default_current_axes
def annotate(artist: Artist, label: str, where: str="over", ax: Axes=None, **kwargs):
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

