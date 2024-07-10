import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
import astropy.units as u
from astropy.time import Time
from scipy.optimize import curve_fit
import all_sky_astro as asa

### step 1 and 3, calculate basic parameters
## f info, location and time, ref catalog, ref image center
f=4.5*(4832/22.3)
time0='2017-05-24 19:00:55'
latitude=38.330278
longitude=74.896667
height=4565
#targetlist='05241900.txt'   # change this for different step
targetlist='ref_XJCAM_20170524190039_1.cat.txt'
x0,y0=np.loadtxt('x0y0.txt',dtype=float)

## read txt
# index / ra / dec / x / y, 5 elements
try:    # for ra dec in float format (degree)
    num,ra1,dec1,x1,y1 = np.loadtxt(targetlist,dtype=float,unpack=True)
    c=SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, frame='icrs')
except: # for ra dec in hmsdms format
    num,ra1,dec1,x1,y1 = np.loadtxt(targetlist,dtype=str,unpack=True)
    c=SkyCoord(ra=ra1, dec=dec1, frame='icrs')
ra=c.ra.degree
dec=c.dec.degree
x=np.array(x1,dtype=float)
y=np.array(y1,dtype=float)

## coordinates
cord = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
loc = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height*u.m)
time = Time(time0)
altaz = cord.transform_to(AltAz(obstime=time,location=loc))

ze=np.pi/2-altaz.alt.rad
az=altaz.az.rad

## initial alt az
# get b',u
b0,u=asa.bu_from_f(x,y,x0,y0,f)
E0,e0=asa.ini_Ee(b0,u,az,ze)

## iteration for basic parameters
for i in range(20):
    E00=E0
    E0,a0,e0=asa.iteration(x,y,x0,y0,az,ze,E0,e0,'iter')
    print(i+1,'E=%.05f, a0=%.05f, e=%.05f, dE:%.03f%%'%(E0,a0,e0,(E00-E0)/E00*100))
Caz,Cze,E,a0,e,k1,k2,k3,k4=asa.iteration(x,y,x0,y0,az,ze,E0,e0,'')
print(k1,k2,k3,k4,E0,a0,e0)
f=open('parameters.txt','w+')
print(E,a0,e,k1,k2,k3,k4,file=f)
f.close()
