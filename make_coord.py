import pytz
import datetime
from timezonefinder import TimezoneFinder
from astropy import units as u
from astropy.coordinates import SkyCoord,AltAz,EarthLocation
from astropy.time import Time
tf = TimezoneFinder()

#for date in YYYYMMDD, time in HH:MM:SS, location in Google Maps format
def return_zenith(date,location,time='12:00:00',frame='icrs'):
    mod_date = f'{date[:4]}-{date[4:6]}-{date[6:]}'
    call_time = Time(f'{mod_date}T{time}')
    loc = EarthLocation.of_address(location)
    geo_loc = loc.to_geodetic()
    tz = tf.timezone_at(lng=geo_loc.lon.value, lat=geo_loc.lat.value)
    timezone = pytz.timezone(tz)
    dt = call_time.to_datetime()
    new_dt = dt - timezone.utcoffset(dt)
    
    coord = SkyCoord(0,90,frame=AltAz(obstime=new_dt.isoformat(),location=loc),unit=u.deg)
    return coord.transform_to(frame)

bday = input("Enter your birthday in YYYYMMDD format:")
bday_time = input("Enter the time you were born in 24-hour format (if unknown, put 12:00)")
bday_loc = input("Enter the location you were born:")

print(return_zenith(bday, bday_loc, bday_time))