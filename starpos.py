from astropy.coordinates import get_moon, AltAz, EarthLocation, get_body
import numpy as np
import astropy.units as u
from astropy.time import Time

def jupiter_position(t):
    height = 2380 * u.m  # m
    longitude = -(16 + 18/60) * u.deg  # deg, negative for west from kyungmin paper
    latitude = (28 + 18/60) * u.deg  # deg from kyungmin paper
    GB_LOCATION = EarthLocation(lat=latitude, lon=longitude, height=height)
    time = Time(t, format='unix', scale='utc')
    jupiter = get_body('jupiter', time, GB_LOCATION).transform_to(AltAz(obstime=time, location=GB_LOCATION))
    az, alt = jupiter.az.value, jupiter.alt.value
    return {"utime":t,"az":az,"alt":alt}

def venus_position(t):
    height = 2380 * u.m  # m
    longitude = -(16 + 18/60) * u.deg  # deg, negative for west from kyungmin paper
    latitude = (28 + 18/60) * u.deg  # deg from kyungmin paper
    GB_LOCATION = EarthLocation(lat=latitude, lon=longitude, height=height)
    time = Time(t, format='unix', scale='utc')
    venus = get_body('venus', time, GB_LOCATION).transform_to(AltAz(obstime=time, location=GB_LOCATION))
    az, alt = venus.az.value, venus.alt.value
    return {"utime":t,"az":az,"alt":alt}



def moon_position(t):
    height = 2380 * u.m  # m
    longitude = -(16 + 18/60) * u.deg  # deg, negative for west from kyungmin paper
    latitude = (28 + 18/60) * u.deg  # deg from kyungmin paper
    GB_LOCATION = EarthLocation(lat=latitude, lon=longitude, height=height)
    time = Time(t, format='unix', scale='utc')
    moon = get_moon(time).transform_to(AltAz(obstime=time, location=GB_LOCATION))
    az, alt = moon.az.value, moon.alt.value
    return {"utime":t,"az":az,"alt":alt}
