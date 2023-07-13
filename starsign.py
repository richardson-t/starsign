from astropy.coordinates import SkyCoord,AltAz,EarthLocation,Angle
from astropy.time import Time
from astropy import units as u
from timezonefinder import TimezoneFinder
from pytz import timezone
from astroquery.simbad import Simbad
from astroquery.hips2fits import hips2fits
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap

class StarSign(object):
    """Information about the star above a given Earth location at a given time.

    Attributes:
        location (str): provided Earth location
        time (astropy.time.Time): provided time, converted to UTC
        coord (astropy.coordinates.SkyCoord): point on the sky over 
            the provided Earth location at the provided time
        star (astropy.table.Table): information about the closest
            star to coord. Retrieved from SIMBAD.
    """
    def __init__(self,location,date,time='12:00:00',frame='icrs'):
        """Instantiate a StarSign.

        Create a StarSign instance given some location, date, time,
        and coordinate frame.
        
        Parameters:
            location (str): Location for the StarSign. Any format 
                that works for Google Maps will function here
                (ex. 'Chicago', 'Chicago, IL', (41.881832, -87.623177)).
            date (str): Day for the StarSign. Should be formatted as 'YYYYMMDD'.
            time (str): Time for the StarSign. Should be formatted using 
                24-hour time as 'HH:MM:SS', or other format that works 
                with astropy.time.Time objects. Defaults to noon.
            frame (str): Reference frame for the StarSign. Can be any of
                astropy's built-in frames ('icrs','fk5','fk4','fk4noeterms','galactic').
                Defaults to 'icrs'.
        """
        self.__location = location
        self.__time = self.__make_time(date,time)
        self.__coord = self.__zenith(location,date,time,frame)
        self.__star = self.__nearest_star()
        
    @property
    def location(self):
        """Returns the StarSign's location attribute."""
        return self.__location
    
    def __interp_date(self,date):
        mod_date = f'{date[:4]}-{date[4:6]}-{date[6:]}'
        return mod_date

    def __make_time(self,date,time):
        mod_date = self.__interp_date(date)
        return Time(f'{mod_date}T{time}')

    @property
    def time(self):
        """Returns the StarSign's time attribute."""
        return self.__time

    def __time_shift(self,dt,eloc):
        geo_loc = eloc.to_geodetic()
        tf = TimezoneFinder()
        tz = tf.timezone_at(lng=geo_loc.lon.value, lat=geo_loc.lat.value)
        zone = timezone(tz)
        return -zone.utcoffset(dt) #negative converts TO UTC, which we want

    def __zenith(self,location,date,time,frame):
        eloc = EarthLocation.of_address(location)
        dt = self.__time.to_datetime()
        dt += self.__time_shift(dt,eloc)
        self.__time = Time(dt.isoformat())
        
        #make SkyCoord directly overhead (90, 0 in alt/az) in location at time
        coord = SkyCoord(0,90,frame=AltAz(obstime=dt.isoformat(),location=eloc),unit=u.deg)
        coord = coord.transform_to(frame)
        return coord

    @property
    def coord(self):
        """Returns the StarSign's coord attribute."""
        return self.__coord

    def __nearest_star(self):
        icrs_coord = self.__coord.transform_to('icrs')
        ra = self.__coord.ra.value
        dec = self.__coord.dec.value
        sign = '+' if dec >= 0 else ''
        stars = Simbad.query_criteria(f'region(circle, ICRS, {ra} {sign}{dec}, 0.5d)',otype='Star')
        return stars[0]

    @property
    def star(self):
        """Returns the StarSign's star attribute."""
        return self.__star

    def visualize(self,hips=None,width=1000,height=1000,
                  fov=0.1*u.deg,coordsys='icrs',cmap='Greys_r'):
        ra = Angle(self.__star['RA'],unit=u.hourangle)
        dec = Angle(self.__star['DEC'],unit=u.deg)
        if hips == None:
            hips = self.__star['MAIN_ID'].split(' ')[0]
        try:
            result = hips2fits.query(hips,width,height,'TAN',
                                     ra,dec,fov,coordsys=coordsys)
        except(AttributeError):
            result = hips2fits.query('Gaia DR3',width,height,'TAN',
                                     ra,dec,fov,coordsys=coordsys)
        im = plt.imshow(result[0].data,cmap=cmap)
        plt.show()
