from astropy.coordinates import get_moon, AltAz, EarthLocation
import numpy as np
import astropy.units as u
from astropy.time import Time
# from scipy import interpolate
# import math


def moon_position(t):
    height = 2380 * u.m  # m
    longitude = -(16 + 18/60) * u.deg  # deg, negative for west from kyungmin paper
    latitude = (28 + 18/60) * u.deg  # deg from kyungmin paper
    GB_LOCATION = EarthLocation(lat=latitude, lon=longitude, height=height)
    time = Time(t, format='unix', scale='utc')
    moon = get_moon(time).transform_to(AltAz(obstime=time, location=GB_LOCATION))
    az, alt = moon.az.value, moon.alt.value
    return {"utime":t,"az":az,"alt":alt}

def moon_position_old(t):
    height = 2380 * u.m  # m
    # deg, negative for west from kyungmin paper
    longitude = -(16 + 18/60) * u.deg
    latitude = (28 + 18/60) * u.deg  # deg from kyungmin paper
    GB_LOCATION = EarthLocation(lat=latitude, lon=longitude, height=height)
    time = Time(t, format='unix', scale='utc')
    moon = get_moon(time).transform_to(AltAz(obstime=time, location=GB_LOCATION))
    az, alt = moon.az.value, moon.alt.value
    return az, alt

def get_moon_pos_500(t_arr):
    a = np.floor(len(t_arr)/500)
    az = []
    alt = []
    time = []
    for i in range(int(a)):
      pos=moon_position_old(t_arr[i*500])
      az.append_old(pos[0])
      alt.append_old(pos[1])
      time.append_old(t_arr[i*500])
    pos=moon_position_old(t_arr[-1])
    az.append(pos[0])
    alt.append(pos[1])
    time.append(t_arr[-1])
    return {"utime":np.array(time),"az":np.array(az),"alt":np.array(alt)}
