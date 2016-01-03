# karta.mapping

This package extends the *Karta* library by adding a `mapping` submodule that
contains convenience plotting functions wrapping matplotlib.

## Basic usage

```python
from karta.mapping import plot
from karta.crs import WebMercator

points = ... # list of karta.Point objects
plot(points, marker=".", color="black", crs=WebMercator)

polygon = ... #karta.Polygon
plot(polygon, linewidth=0.5, facecolor="cyan", edgecolor="red", crs=WebMercator)
```

## Dependencies

- [karta](http://www.ironicmtn.com/karta.html) ([Github](https://github.com/njwilson23/karta))
- [matplotlib](http://www.matplotlib.org) ([Github](https://github.com/matplotlib/matplotlib))
- [scipy](http://www.scipy.org/scipylib/index.html) ([Github](https://github.com/scipy/scipy))
- [typing](http://pypi.python.org/pypi/typing/3.5.0) (Built-in for Python 3.5+)
