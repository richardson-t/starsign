import pytest
from astropy.time import Time
from astropy.coordinates import Angle
from astropy import units as u
from timezonefinder import TimezoneFinder
from pytz import timezone
from starsign.starsign import StarSign

locations = ['Chicago, IL', 'Albany, USA']
dates = ['2023-07-13','2023-01-13','2020-02-29']
times = ['20:00', '00:00']
frames = ['icrs','fk5','fk4','fk4noeterms','galactic']

#unit tests

def test_functioning():
    for loc in locations:
        s = StarSign(loc,dates[0])
    for date in dates:
        s = StarSign(locations[0],date)
    for time in times:
        s = StarSign(locations[0],dates[0],time=time)
    for frame in frames:
        s = StarSign(locations[0],dates[0],frame=frame)

def test_time_creation():
    #does it read in and shift time appropriately
    for i in range(2):
        t = Time(f'{dates[i]}T{times[i]}')
        dt = t.to_datetime()
        tf = TimezoneFinder()
        tz = tf.timezone_at(lng=-87.6298, lat=41.8781)
        zone = timezone(tz)
        offset = zone.utcoffset(dt)
        
        s = StarSign(locations[0],dates[i],times[i])
        sign_delta = -((s.time-t)*86400).value

        assert offset.total_seconds() == pytest.approx(sign_delta,abs=1e-3)
    
def test_frame():
    #are the coordinates in the right frames? are they the right numbers?
    testers = ['nonsense',None,0.5]
    for test in testers:
        try:
            s = StarSign(locations[0],dates[0],frame=test)
        except(ValueError):
            continue
    
def test_star_finding():
    #does it return the same star for different objects with the same input?
    #is the star close to the expected coordinates?
    s1 = StarSign(locations[0],dates[0])
    s2 = StarSign(locations[0],dates[0])
    assert s1.star['MAIN_ID'] == s2.star['MAIN_ID']

    ra = Angle(s1.star['RA'],unit=u.hourangle).to(u.deg)
    dec = Angle(s1.star['DEC'],unit=u.deg)

    assert ra.value == pytest.approx(s1.coord.ra.value,abs=0.25)
    assert dec.value == pytest.approx(s1.coord.dec.value,abs=0.25)
    
#end-to-end tests

def test_vis():
    #does it produce a plot that is consistent with input?
    #note: this requires by-eye validation
    s = StarSign(locations[0],dates[0])
    s.visualize()
    s.visualize(hips='Gaia DR3')
    s.visualize(width=3000)
    s.visualize(height=3000)
    s.visualize(fov=1)
    s.visualize(cmap='viridis')

    #does it handle edge cases?
    #image too big to get from HiPS
    try:
        s.visualize(width=10000,height=10000)
    except(RuntimeError):
        pass

    #one of the dimensions is zero
    try:
        s.visualize(width=0)
    except(ValueError):
        pass
    try:
        s.visualize(height=0)
    except(ValueError):
        pass

    
