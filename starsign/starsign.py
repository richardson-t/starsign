from astropy.coordinates import SkyCoord,AltAz,EarthLocation
from astropy.time import Time
from astropy import units as u
from timezonefinder import TimezoneFinder
from pytz import timezone
from astroquery.simbad import Simbad

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
        stars = Simbad.query_criteria(f'region(circle, ICRS, {ra} {sign}{dec}, 0.5d)',maintype='Star', otypes=['Ma*','MS*','Ev*','**','LM*','V*','Em*','PM*','HV*'])
        return stars[0]

    @property
    def star(self,radius='0d20m0s'):
        return self.__star

otype_d = {
    'Ma*' : 'Massive Star',
    'bC*' : 'Beta Cepheid Variable',
    'sg*' : 'Evolved Supergiant',
    's*r' : 'Red Supergiant',
    's*y' : 'Yellow Supergiant',
    's*b' : 'Blue Supergiant',
    'WR*' : 'Wolf-Rayet',
    'N*' : 'Neutron Star',
    'Psr' : 'Pulsar',
    'MS*' : 'Main Sequence Star',
    'Be*' : 'Be Star',
    'BS*' : 'Blue Straggler',
    'gD*' : 'Gamma Doradus Variable Star',
    'dS*' : 'Delta Scuti Variable Star',
    'Ev*' : 'Evolved Star',
    'RG*' : 'Red Giant Branch Star',
    'HS*' : 'Hot Subdwarf',
    'HB*' : 'Horizontal Branch Star',
    'RR*' : 'RR Lyrae Variable Star',
    'WV*' : 'Type II Cepheid Variable Star',
    'Ce*' : 'Cepheid Variable Star',
    'cC*' : 'Classical Cepheid Variable Star',
    'C*' : 'Carbon Star',
    'S*' : 'S-type Star',
    'LP*' : 'Long-Period Variable Star',
    'AB*' : 'Asymptotic Giant Branch Star',
    'Mi*' : 'Mira Variable Star',
    'OH*' : 'OH/IR Star',
    'pA*' : 'Post-AGB Star',
    'RV*' : 'RV Tauri Variable Star',
    'PN' : 'Planetary Nebula',
    'WD*' : 'White Dwarf',
    '**' : 'Double or Multiple Star',
    'El*' : 'Ellipsoidal Variable',
    'EB*' : 'Eclipsing Binary',
    'SB*' : 'Spectroscopic Binary',
    'BY*' : 'BY Draconis Variable',
    'RS*' : 'RS Canum Venaticorum Variable',
    'Sy*' : 'Symbiotic Star',
    'XB*' : 'X-ray Binary',
    'LXB' : 'Low Mass X-ray Binary',
    'HXB' : 'High Mass X-ray Binary',
    'CV*' : 'Cataclysmic Binary',
    'No*' : 'Classical Nova',
    'SN*' : 'Supernova',
    'LM*' : 'Low Mass Star',
    'BD*' : 'Brown Dwarf',
    'V*' : 'Variable Star',
    'Ir*' : 'Irregular Variable Star',
    'Er*' : 'Eruptive Variable Star',
    'Ro*' : 'Rotating Variable Star',
    'Pu*' : 'Pulsating Variable Star',
    'Em*' : 'Emission-line Star',
    'PM*' : 'High Proper Motion Star',
    'HV*' : 'High Velocity Star'
}