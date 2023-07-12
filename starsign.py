from astropy.coordinates import SkyCoord,AltAz,EarthLocation
from astropy.time import Time
from astropy import units as u
from timezonefinder import TimezoneFinder
from pytz import timezone
from astroquery.simbad import Simbad

class StarSign(object):    
    def __init__(self,location,date,time='12:00:00',frame='icrs'):
        self.location = location
        self.date = date
        self.time = time
        self.coord = self._zenith(location,date,time,frame)
        
    def _interp_date(self,date):
        mod_date = f'{date[:4]}-{date[4:6]}-{date[6:]}'
        return mod_date

    def _make_time(self,date,time):
        mod_date = self._interp_date(date)
        return Time(f'{mod_date}T{time}')
        
    def _zenith(self,location,date,time,frame):
        #read in provided location and time
        loc = EarthLocation.of_address(location)
        dt = self._make_time(date,time).to_datetime()

        #use location to convert to UTC
        geo_loc = loc.to_geodetic()
        tf = TimezoneFinder()
        tz = tf.timezone_at(lng=geo_loc.lon.value, lat=geo_loc.lat.value)
        zone = timezone(tz)
        dt -= zone.utcoffset(dt) #utcoffset works FROM UTC, but we want TO
        
        #make SkyCoord directly overhead (90, 0 in alt/az) in location at time
        coord = SkyCoord(0,90,frame=AltAz(obstime=dt.isoformat(),location=loc),unit=u.deg)
        coord = coord.transform_to(frame)
        return coord

    def object(self,radius='0d20m0s'):
        return Simbad.query_region(self.coord,radius)[0]
