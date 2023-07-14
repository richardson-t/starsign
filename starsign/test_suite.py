import pytest
from starsign import StarSign

locations = ['Chicago, IL', 'Albany, USA']
dates = ['2023-07-13', '2020-02-29']
times = ['20:00', '00:00']
frames = ['icrs','fk5','fk4','fk4noeterms','galactic']

#unit tests

def test_time_creation():
    #does it read in and shift time appropriately

def test_coordinate_creation():
    #are the coordinates in the right frames? are they the right numbers?
    
def test_star_finding():
    #does it return the same star for different objects with the same input?
    #is the star close to the expected coordinates?

#end-to-end tests

def test_vis():
    #does it produce a plot that is consistent with input?
    #note: this requires by-eye validation
    s = StarSign(locations[0],dates[0])
    s.visualize()
    s.visualize(hips='Gaia DR3')
    s.visualize(width=3000)
    s.visualize(height=3000)
    s.visualize(fov=10)
    s.visualize(cmap='viridis')

    #does it handle edge cases?
    #image too big to get from HiPS
    try:
        s.vizualize(width=10000,height=10000)
    except(RuntimeError):
        continue
    
