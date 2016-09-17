"""
mapping is a module for plotting vector and raster objects from the karta
package using matplotlib
"""

import functools

import scipy.optimize
import numpy as np

from typing import Union, Iterable, Tuple, Callable
from matplotlib.pyplot import gca, sci, Axes, Artist, cm

from karta.vector import Point, Line, Polygon, Geometry
from karta.vector import Multipoint, Multiline, Multipolygon
from karta.raster import RegularGrid
from karta.crs import CRS, Cartesian, SphericalEarth

def default_current_axes(wrappedfunc: Callable):
    """ Decorator to set current Axes as default ax in plotting functions """
    @functools.wraps(wrappedfunc)
    def replacementfunc(*args, **kwargs):
        if "ax" not in kwargs:
            # don't use dict.setdefault because don't want gca() to be
            # evaluated unconditionally
            kwargs["ax"] = gca()
        return wrappedfunc(*args, **kwargs)
    return replacementfunc

def recurse_iterables(wrappedfunc: Callable):
    """ Decorator to generate functions that apply themselve recursively to
    iterable non-Geometry inputs. """
    @functools.wraps(wrappedfunc)
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

def get_axes_limits(ax: Axes, ax_crs: CRS, crs: CRS=SphericalEarth):
    """ Get the limits of the window covered by an Axes in another coordinate
    system. """
    xl, xr = ax.get_xlim()
    yb, yt = ax.get_ylim()
    ax_bbox = get_axes_extent(ax, ax_crs, crs=crs)

    # Minimize bottom spine
    x_ = scipy.optimize.fminbound(lambda x: ax_crs.transform(crs, x, yb)[1], xl, xr)
    ymin = ax_crs.transform(crs, x_, yb)[1]

    # Maximize top spine
    x_ = scipy.optimize.fminbound(lambda x: -ax_crs.transform(crs, x, yt)[1], xl, xr)
    ymax = ax_crs.transform(crs, x_, yt)[1]

    # Minimize left spine
    y_ = scipy.optimize.fminbound(lambda y: ax_crs.transform(crs, xl, y)[0], yb, yt)
    xmin = ax_crs.transform(crs, xl, y_)[0]

    # Maximize right spine
    y_ = scipy.optimize.fminbound(lambda y: -ax_crs.transform(crs, xr, y)[0], yb, yt)
    xmax = ax_crs.transform(crs, xr, y_)[0]
    return xmin, xmax, ymin, ymax

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
def add_graticule(xs: Iterable[float], ys: Iterable[float], ax: Axes=None,
        map_crs: CRS=Cartesian, graticule_crs: CRS=SphericalEarth,
        nsegments=25, lineargs=None):
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
    nsegments : int, optional
        Number of segments to use to approximate curving graticule lines
    lineargs : dict, optional
        Arguments passed to karta.mapping.plot while drawing graticule lines
    """
    if lineargs is None:
        lineargs = dict()

    lineargs.setdefault("color", "black")
    lineargs.setdefault("linewidth", 0.5)

    xmin, xmax, ymin, ymax = get_axes_limits(ax, map_crs, crs=graticule_crs)
    artists = []
    for i in range(len(xs)):
        coords = zip(np.linspace(xs[i], xs[i], nsegments+1),
                     np.linspace(ymin, ymax, nsegments+1))
        _line = Line(coords, crs=graticule_crs)
        artists.append(plot(_line, ax=ax, crs=map_crs, **lineargs))
    for i in range(len(ys)):
        coords = zip(np.linspace(xmin, xmax, nsegments+1),
                     np.linspace(ys[i], ys[i], nsegments+1))
        _line = Line(coords, crs=graticule_crs)
        artists.append(plot(_line, ax=ax, crs=map_crs, **lineargs))
    return artists

def isbetween(x: float, a: float, b: float) -> bool:
    return (a < x < b) or (b < x < a)

def froot(f: float, a: float, b: float) -> float:
    return scipy.optimize.brentq(f, a, b)

@default_current_axes
def label_ticks(xs: Iterable[float], ys: Iterable[float], ax: Axes=None,
        map_crs: CRS=Cartesian, graticule_crs: CRS=SphericalEarth,
        textargs=None, tickargs=None,
        xformatter=None, yformatter=None):
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
    xformatter : callable, optional
        function that given an easting/longitude returns a label
    yformatter : callable, optional
        function that given a northing/latitude returns a label
    """
    if textargs is None:
        textargs = dict()

    if tickargs is None:
        tickargs = dict(marker="+", mew=2, ms=14, mfc="k", mec="k", ls="none")

    if xformatter is None:
        xformatter = lambda x: "{0} E".format(x)

    if yformatter is None:
        yformatter = lambda y: "{0} N".format(y)

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
                                          xmin, xmax), ymin, xformatter(x)))

    for y in ys:
        if isbetween(y, bbox[0][1], bbox[1][1]):
            ticks["yticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymin)[1]-y,
                                          xmin, xmax), ymin, yformatter(y)))

    # top spine
    for x in xs:
        if isbetween(x, bbox[2][0], bbox[3][0]):
            ticks["xticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymax)[0]-x,
                                          xmin, xmax), ymax, xformatter(x)))

    for y in ys:
        if isbetween(y, bbox[2][1], bbox[3][1]):
            ticks["yticks"].append((froot(lambda xt:
                                          map_crs.transform(graticule_crs, xt, ymax)[1]-y,
                                          xmin, xmax), ymax, yformatter(y)))

    # left spine
    for x in xs:
        if isbetween(x, bbox[0][0], bbox[3][0]):
            ticks["xticks"].append((xmin,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmin, yt)[0]-x,
                                          ymin, ymax), xformatter(x)))


    for y in ys:
        if isbetween(y, bbox[0][1], bbox[3][1]):
            ticks["yticks"].append((xmin,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmin, yt)[1]-y,
                                          ymin, ymax), yformatter(y)))


    # right spine
    for x in xs:
        if isbetween(x, bbox[1][0], bbox[2][0]):
            ticks["xticks"].append((xmax,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmax, yt)[0]-x,
                                          ymin, ymax), xformatter(x)))

    for y in ys:
        if isbetween(y, bbox[1][1], bbox[2][1]):
            ticks["yticks"].append((xmax,
                                    froot(lambda yt:
                                          map_crs.transform(graticule_crs, xmax, yt)[1]-y,
                                          ymin, ymax), yformatter(y)))

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

def plot(geom: Union[Union[Geometry, RegularGrid], Iterable[Union[Geometry, RegularGrid]]], *args, **kwargs):
    """ Metafunction that dispatches to the correct plotting routine. """
    if isinstance(geom, list):
        results = []
        for g in geom:
            results.append(plot(g, *args, **kwargs))
        return results
    else:
        if isinstance(geom, RegularGrid):
            func = plot_grid
        elif not hasattr(geom, "_geotype"):
            raise TypeError("Invalid input type: {0}".format(type(geom)))
        elif geom._geotype == "Point":
            func = plot_point
        elif geom._geotype == "Line":
            func = plot_line
        elif geom._geotype == "Polygon":
            func = plot_polygon
        elif geom._geotype == "Multipoint":
            func = plot_multipoint
        elif geom._geotype == "Multiline":
            func = plot_multiline
        elif geom._geotype == "Multipolygon":
            func = plot_multipolygon
        else:
            raise TypeError("Invalid geotype: {0}".format(geom._geotype))
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
def plot_multipoint(geom: Union[Multipoint, Iterable[Multipoint]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("linestyle", "none")
    kwargs.setdefault("marker", ".")
    x, y = geom.get_coordinate_lists(crs=crs)
    return ax.plot(x, y, *args, **kwargs)

@default_current_axes
def plot_multiline(geom: Union[Multiline, Iterable[Multiline]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    out = []
    for line in geom:
        x, y = line.get_coordinate_lists(crs=crs)
        out.append(ax.plot(x, y, *args, **kwargs))
    return out

@default_current_axes
def plot_multipolygon(geom: Union[Multipolygon, Iterable[Multipolygon]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("facecolor", "none")
    kwargs.setdefault("edgecolor", "black")
    out = []
    for polygon in geom:
        x, y = polygon.get_coordinate_lists(crs=crs)
        out.append(ax.fill(x, y, *args, **kwargs))
    return out

@default_current_axes
def plot_grid(grid: RegularGrid, ax: Axes=None, crs: CRS=None, band: Union[int, tuple]=-1, **kwargs):
    """ Plot a grid instance

    Parameters
    ----------
    grid : RegularGrid
        raster data to plot
    ax : Axes, optional
        Axes to plot to [default plt.gca()]
    crs : CRS, optional
        Currently not supported
    band : int or tuple, optional
        Band(s) to plot. If *grid* has three bands, by default the three are
        plotted in false colour as RGB channels. Otherwise, the first band is
        plotted by default. If *band* is a tuple, it must have three integer
        elements.

    Notes
    -----
    Additional arguments are passed to `matplotlib.pyplot.imshow`
    """
    kwargs.setdefault("origin", "bottom")
    kwargs.setdefault("extent", grid.get_extent(crs=crs))
    kwargs.setdefault("cmap", cm.binary_r)

    if crs is not None and crs != grid.crs:
        raise NotImplementedError("RegularGrid reprojection not supported")

    # compute the pixels that can actually be displayed
    # be slightly generous by using a factor of 0.75 to avoid choosing too low
    # of a resolution
    _, _, width, height = ax.bbox.bounds
    ny, nx = grid.size
    r = (max(int(0.75*ny//height), 1), max(int(0.75*nx//width), 1))
    if band == -1:
        if len(grid.bands) == 3 and (band == -1):
            band = (0, 1, 2)
        else:
            band = 0
    if isinstance(band, int):
        arr = grid[::r[0],::r[1],band]
        arr = np.ma.masked_equal(arr, grid.nodata)
    else:
        if len(band) not in (3, 4):
            raise ValueError("bands must be RGB or RGBA (length 3 or 4)")
        arr = np.dstack([grid[::r[0],::r[1],i] for i in band]).astype(np.float32)
        arr = np.ma.masked_equal(arr, grid.nodata)
        arr[:,:,:3] /= arr[:,:,:3].max()

    im = ax.imshow(arr, **kwargs)
    if ax == gca():
        sci(im)
    return im

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

