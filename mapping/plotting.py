import functools
from typing import Union, Iterable, Tuple, Callable

import numpy as np
from karta.vector import Point, Line, Polygon, Geometry
from karta.vector import Multipoint, Multiline, Multipolygon
from karta.raster import RegularGrid
from karta.crs import CRS

import matplotlib.path
import matplotlib.collections
from matplotlib.pyplot import gca, sci, Axes, Artist, cm

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
    kwargs.setdefault("facecolors", "none")
    kwargs.setdefault("edgecolors", "black")
    paths = [matplotlib.path.Path(vertices.asarray()[:,:2], readonly=True)
            for vertices in geom.get_vertices(crs=crs)]
    coll = matplotlib.collections.PathCollection(paths, *args, **kwargs)
    ax.add_artist(coll)
    return coll

@default_current_axes
def plot_multipolygon(geom: Union[Multipolygon, Iterable[Multipolygon]], *args,
        ax: Axes=None, crs: CRS=None, **kwargs):
    """ Plot a Line geometry, projected to the coordinate system `crs` """
    kwargs.setdefault("facecolors", "none")
    kwargs.setdefault("edgecolors", "black")
    paths = [matplotlib.path.Path(vertices.asarray()[:,:2], closed=True, readonly=True)
            for vertices in geom.get_vertices(crs=crs)]
    coll = matplotlib.collections.PathCollection(paths, *args, **kwargs)
    ax.add_artist(coll)
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

