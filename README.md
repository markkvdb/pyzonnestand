# PyZonnestand

Fast and accurate calculation of the position of the sun.

## Quickstart

```python
from pyzonnestand import sun_pos

start = datetime(2019, 1, 1, 12, 0, 0)
dt = [start + timedelta(hours=i) for i in range(48)]
data = sunpos(
    dt=dt,
    latitude=52.0,
    longitude=5.0,
    elevation=0.0,
)
```


## Credits

Heavily inspired by [sun-position](https://github.com/s-bear/sun-position).
