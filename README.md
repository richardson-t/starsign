# starsign

Given an Earth location, date, and time, this package can find the star closest to directly overhead. If you want a horoscope backed by the full might of SIMBAD, then it's your lucky day!

**Example usage:**
```
from starsign.starsign import StarSign
s = StarSign('Chicago, IL','2023-07-14')
```
StarSign will transform the provided location on Earth into a coordinate in the sky and find the closest object identified as a star in [SIMBAD](https://simbad.u-strasbg.fr/simbad/sim-fbasic). Once this is done, SIMBAD's information on the object can be accessed:
```
print(s.coord)
<SkyCoord (ICRS): (ra, dec) in deg
    (99.27076837, 41.89516319)>
print(s.star['MAIN_ID'])
BD+41  1471
print(s.star['OTYPE'])
HighPM*
```
The object and surrounding area can be visualized by running:
```
s.visualize()
```
For more details on functionality, run help(StarSign) from within a Python session. This software package can be installed by running:
```
pip install starsign
```
from a terminal window or downloaded from [our PyPI page](https://pypi.org/project/starsign/).
