from astropy.coordinates import SkyCoord,AltAz,EarthLocation
from astropy.time import Time
from astropy import units as u
from timezonefinder import TimezoneFinder
from pytz import timezone
from astroquery.simbad import Simbad
custom_simbad = Simbad()
custom_simbad.add_votable_fields('otype','otype(opt)', 'otypes')

class StarSign(object):    
    def __init__(self,location,date,time='12:00:00',frame='icrs'):
        self.__location = location
        self.__time = self.__make_time(date,time)
        self.__coord = self.__zenith(location,date,time,frame)
        self.__star = self.__nearest_star()
        
    @property
    def location(self):
        return self.__location
    
    def __interp_date(self,date):
        mod_date = f'{date[:4]}-{date[4:6]}-{date[6:]}'
        return mod_date

    def __make_time(self,date,time):
        mod_date = self.__interp_date(date)
        return Time(f'{mod_date}T{time}')

    @property
    def time(self):
        return self.__time

    def __time_shift(self,dt,eloc):
        geo_loc = eloc.to_geodetic()
        tf = TimezoneFinder()
        tz = tf.timezone_at(lng=geo_loc.lon.value, lat=geo_loc.lat.value)
        zone = timezone(tz)
        return -zone.utcoffset(dt) #negative converts TO UTC, which we want

    def __zenith(self,location,date,time,frame):
        #read in provided location and time
        eloc = EarthLocation.of_address(location)
        dt = self.__time.to_datetime()
        dt += self.__time_shift(dt,eloc)
        #make SkyCoord directly overhead (90, 0 in alt/az) in location at time
        coord = SkyCoord(0,90,frame=AltAz(obstime=dt.isoformat(),location=eloc),unit=u.deg)
        coord = coord.transform_to(frame)
        return coord

    @property
    def coord(self):
        return self.__coord

    def __nearest_star(self):
        ra = self.__coord.ra.value
        dec = self.__coord.dec.value
        sign = '+' if dec >= 0 else ''
        stars = custom_simbad.query_criteria(f'region(circle, ICRS, {ra} {sign}{dec}, 0.5d)', otype='Star')
        return stars[0]

    @property
    def star(self,radius='0d20m0s'):
        return self.__star
