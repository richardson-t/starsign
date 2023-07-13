from astropy.coordinates import SkyCoord,AltAz,EarthLocation,Angle
from astropy.time import Time
from astropy import units as u
from timezonefinder import TimezoneFinder
from pytz import timezone
from astroquery.simbad import Simbad
from astroquery.hips2fits import hips2fits
import numpy as np
import matplotlib.pyplot as plt

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
        self.__star_info = self.__nearest_star()
        
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
        
        custom_simbad = Simbad()
        custom_simbad.add_votable_fields('otype')
        stars = custom_simbad.query_criteria(f'region(circle, ICRS, {ra} {sign}{dec}, 0.5d)', otype='Star')
        return stars[0]

    @property
    def star(self):
        """Returns the StarSign's star attribute."""
        return self.__star

    def __width_long(self,width,height):
        if width >= height:
            return True
        else:
            return False
    
    def visualize(self,hips=None,width=1000,height=1000,
                  fov=0.1*u.deg,cmap='Greys_r'):
        ra = Angle(self.__star['RA'],unit=u.hourangle)
        dec = Angle(self.__star['DEC'],unit=u.deg)
        if hips == None:
            hips = self.__star['MAIN_ID'].split(' ')[0]

        try:
            #sometimes object identifiers don't have associated catalogs;
            #if not, default to Gaia
            try:
                result = hips2fits.query(hips,width,height,'TAN',
                                         ra,dec,fov)
            except(AttributeError):
                result = hips2fits.query('Gaia DR3',width,height,'TAN',
                                         ra,dec,fov)
        except:
            raise RuntimeError('Read from HiPS timed out; try fewer pixels')
                
        im = plt.imshow(result[0].data,cmap=cmap)

        min_px,max_px = min(width,height),max(width,height)
        ticks = np.linspace(.25,.75,3)*max_px
        short_ticks = np.linspace(min(ticks),max(ticks),len(ticks))-(max_px-min_px)/2
        short_ticks = short_ticks[np.logical_and(short_ticks > 0,short_ticks < min_px)]
        tick_scale = np.format_float_positional(fov.value/4,precision=3)
        pm = f'{tick_scale}'
        
        if self.__width_long(width,height):
            xlabels = [f'-{pm}',f"{np.round(ra.to(u.deg).value,3)}",f'+{pm}']
            ylabels = [None,f'{np.round(dec.value,3)}',None]
            if len(short_ticks) < len(ticks):
                ylabels = [ylabels[1]]
            plt.xticks(ticks=ticks,labels=xlabels)
            plt.yticks(ticks=short_ticks,labels=ylabels)
        else:
            xlabels = [None,f'{np.round(ra.to(u.deg).value,3)}',None]
            ylabels = [f'-{pm}',f'{np.round(dec.value,3)}',f'+{pm}']
            if len(short_ticks) < len(ticks):
                xlabels = [xlabels[1]]
            plt.xticks(ticks=short_ticks,labels=xlabels)
            plt.yticks(ticks=ticks,labels=ylabels)

        plt.gca().tick_params(width=2)
        plt.title(f"{self.__star['MAIN_ID']}")
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.show()
