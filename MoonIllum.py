from __future__ import print_function, division
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_moon, SkyCoord
from astropy.time import Time
from datetime import date, datetime, time,timedelta
import numpy as np
from astropy.table import Table
import time
from PyAstronomy import pyasl 
import ephem
import astroplan
import requests
from requests.auth import HTTPBasicAuth
import astropy.units as u
import pytz
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import astropy.units as u
from astropy.coordinates import *
import math
from astropy import coordinates as coord
import skyfield
from skyfield.api import Topos, load, Star
from skyfield.toposlib import Topos
import sys

names=['2019cwt']

tomorrow = datetime.now() + timedelta(days = (1))
timedate = tomorrow.replace(hour=4)
timedate = timedate.replace(minute=0)
timedate = timedate.replace(second=0)
print(timedate)


ts = load.timescale()
DEG_IN_RADIAN = 57.2957795130823
EARTHRAD_IN_AU = 23454.7910556298
HRS_IN_RADIAN = 3.81971863420549

def query_mag(name):
    url = 'https://ziggy.ucolick.org/yse/download_photometry/' + name + '/'
    rr = requests.get(url=url, auth=HTTPBasicAuth('djones','BossTent1'))
    fout = open(name+'tmp','w')
    print(rr._content.decode('utf-8'),file=fout)
MIN_IN_DEG = SEC_IN_MIN = 60
def ephem_angle_to_float(angle):
    degrees, minutes, seconds = [float(s) for s in str(angle).split(':')]
    value = abs(degrees) + \
        minutes / MIN_IN_DEG + \
        seconds / SEC_IN_MIN / MIN_IN_DEG
    return math.copysign(value, degrees)
def JulianDay (date,month,year,UT):
	if (month<=2):
		month=month+12; year=year-1
	return int(365.25*year) + int(30.6001*(month+1)) - 15 + 1720996.5 + date + UT/24.0
def lunskybright(alpha,rho,kzen,altmoon,alt,moondist):
	rho_rad = rho/DEG_IN_RADIAN
	alpha = (180. - alpha)
	Zmoon = (90. - altmoon)/DEG_IN_RADIAN
	Z = (90. - alt)/DEG_IN_RADIAN
	moondist = EARTHRAD_IN_AU * moondist/(60.27)
	istar = -0.4*(3.84 + 0.026*np.abs(alpha) + 4.0e-9*pow(alpha,4.))
	istar =  pow(10.,istar)/(moondist * moondist)
	if(np.abs(alpha) < 7.):
		istar = istar * (1.35 - 0.05 * np.abs(istar))
	fofrho = 229087. * (1.06 + np.cos(rho_rad)*np.cos(rho_rad))
	if(np.abs(rho) > 10.):
		fofrho=fofrho + pow(10.,(6.15 - rho/40.))
	elif (np.abs(rho) > 0.25):
		fofrho= fofrho+ 6.2e7 / (rho*rho)
	else:
		fofrho = fofrho+9.9e8
	Xzm = np.sqrt(1.0 - 0.96*np.sin(Zmoon)*np.sin(Zmoon))
	if(Xzm != 0.):
		Xzm = 1./Xzm
	else:
		Xzm = 10000.
	Xo = np.sqrt(1.0 - 0.96*np.sin(Z)*np.sin(Z))
	if(Xo != 0.):
		Xo = 1./Xo
	else:
		Xo = 10000.
	Bmoon = fofrho * istar * pow(10.,(-0.4*kzen*Xzm)) * (1. - pow(10.,(-0.4*kzen*Xo)))
	#f(Bmoon > 0.001):
	return(22.50 - 1.08574 * np.log(Bmoon/34.08))
	#else:
	#	return(99.)
def decdeg2dms(dd):
   is_positive = dd >= 0
   dd = abs(dd)
   minutes,seconds = divmod(dd*3600,60)
   degrees,minutes = divmod(minutes,60)
   degrees = degrees if is_positive else -degrees
   return (degrees,minutes,seconds)
eph = load('de421.bsp')
moon = eph['moon']
sun = eph['sun']
earth = eph['earth']

swope = earth + Topos('29.01597 S','70.69208 W')
moonobs = swope.at(ts.utc(timedate.year, timedate.month,timedate.day,timedate.hour,0,0)).observe(moon)
sunobs = swope.at(ts.utc(timedate.year, timedate.month,timedate.day,timedate.hour,0,0)).observe(sun)
moonapp = moonobs.apparent()
sunapp = sunobs.apparent()

elongation = (moonapp.separation_from(sunapp)).to(u.deg)

print("Angle from Moon to Sun: " + str(elongation))
###########################################################################

alt, az, distance = moonapp.altaz()
moonapp.radec()
ra, dec, distance = moonapp.radec()
print(ra.hstr())
print(dec.dstr())
print("Altitude of the moon: " + str(alt))
print("Distance to the moon: " + str(distance))
###########################################################################
#Finding the altitude of the object
for name in names:
	print("STARTING " + name)
	query_mag(name)
	txt_file = (name+"tmp")
	old="s"
	RALine=""
	DecLine=""

	with open(name+'tmp', 'r') as f:
	    content = f.readlines()
	    for line in content:
	    	if line.find("# RA:") != -1:
	    		RALine=line
	    	if line.find("# DEC:") != -1:
	    		DecLine=line
	if DecLine == "" or RALine =="":
		print("Issue with the tmp file.")

	DECobject=float(DecLine[7:])
	RAobject=float(RALine[6:])
	print("RA, DEC of object: " + str(RAobject) + ',' + str(DECobject))

	RAhrs=decdeg2dms(RAobject)
	DEChrs = decdeg2dms(DECobject)

	target = Star(ra_hours=(RAobject/360 * 24,0,0),dec_degrees=(DEChrs[0],DEChrs[1],DEChrs[2]))
	objecto = swope.at(ts.utc(timedate.year, timedate.month,timedate.day,timedate.hour,0,0)).observe(target)
	alto, azo, do = objecto.apparent().altaz()

	print("Altitude of the object: " + str(alto))
	###########################################################################

	rho = (moonapp.separation_from(objecto.apparent())).to(u.deg)
	print("angle of moon and object: " + str(rho))

	kzen = 0.172

	elong = str(elongation.to(u.rad).to(u.deg))
	altMoon = str(alt.to(u.deg))
	alt_target = str(alto.to(u.deg))
	dist = str(distance)
	brightness = lunskybright(float(elong[:elong.index(' ')]), float(str(rho)[:str(rho).index(' ')]), kzen, float(altMoon[:altMoon.index(' ')]), float(alt_target[:alt_target.index(' ')]), float(dist[:dist.index(' ')]))
	print("Lunar Sky Brightness of " + name + " is " + str(brightness))
	print("ENDING " + name)