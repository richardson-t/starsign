from astropy import units as u
from astropy.coordinates import SkyCoord,AltAz,EarthLocation
from astropy.time import Time

#for date in YYYYMMDD, time in HH:MM:SS, location in Google Maps format
def return_zenith(date,location,time='12:00:00',frame='icrs'):
    mod_date = f'{date[:4]}-{date[4:6]}-{date[6:]}'
    call_time = Time(f'{mod_date}T{time}')
    loc = EarthLocation(location)
    
    coord = SkyCoord(0,90,frame=AltAz(obstime=call_time,location=loc),unit=u.deg)
    return coord.transform_to(frame)
