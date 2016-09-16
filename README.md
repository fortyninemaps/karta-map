# mapping

This module adds [*matplotlib*](https://github.com/matplotlib/matplotlib)
plotting functions to the [*Karta* geospatial library](https://github.com/fortyninemaps/karta).

## Basic usage

```python
from mapping import plot, annotate
from karta.crs import WebMercator

points = ...    # list of karta.Point objects
plot(points, marker=".", color="black", crs=WebMercator)

polygon = ...   # karta.Polygon
artist = plot(polygon, linewidth=0.5, facecolor="cyan", edgecolor="red", crs=WebMercator)
annotate(artist[0], "Study region", where="over")
```

## Dependencies

`mapping` is written in type-annotated Python, and therefore requires Python 3.

- [karta](http://www.fortyninemaps.com/karta.html) ([Github](https://github.com/njwilson23/karta))
- [matplotlib](http://www.matplotlib.org) ([Github](https://github.com/matplotlib/matplotlib))
- [scipy](http://www.scipy.org/scipylib/index.html) ([Github](https://github.com/scipy/scipy))

Versions of Python<3.5 additionally require the [typing](https://pypi.python.org/pypi/typing/) module.
