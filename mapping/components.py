from typing import Iterable

import numpy as np
import scipy.optimize
from karta.vector import Point, Line, Polygon
from karta.crs import CRS, Cartesian, SphericalEarth

from matplotlib.pyplot import Axes

from .plotting import plot, default_current_axes

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

