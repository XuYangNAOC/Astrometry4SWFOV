from astropy.io import fits
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from scipy.optimize import curve_fit
import all_sky_astro as asa

### step 3, get parameters
## get basic info
# parameters
x0,y0=np.loadtxt('x0y0.txt',dtype=float)
f=4.5*(4832/22.3)

# time
time0='2017-05-24 19:00:55'
time = Time(time0)

# location
latitude=38.330278
longitude=74.896667
height=4565
loc = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height*u.m)

# target lists
targetlist='ref_XJCAM_20170524190039_1.cat.txt'

# loop length and step
length=20
step=0.1

## read cross-matched catalog
# target catalog name
num,ra1,dec1,x1,y1 = np.loadtxt(targetlist,dtype=float,unpack=True)
c=SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree, frame='icrs')
ra=c.ra.degree
dec=c.dec.degree
x=np.array(x1,dtype=float)
y=np.array(y1,dtype=float)

## fit new x0y0
length=20
step1=1
step2=0.1
X,Y,A,C=asa.findcenter(ra,dec,x,y,x0,y0,time0,loc,f,length,step1)
#X,Y,A,C=np.loadtxt('findcenter.txt',dtype=float,unpack=True)
x0_new1,y0_new1=asa.fitcenter(X,Y,A,C,length,step1)

X,Y,A,C=asa.findcenter(ra,dec,x,y,x0_new1,y0_new1,time0,loc,f,length,step2)
#X,Y,A,C=np.loadtxt('findcenter.txt',dtype=float,unpack=True)
x0_new,y0_new=asa.fitcenter(X,Y,A,C,length,step2)

# save
with open('x0y0_new.txt','w+') as files:
    print(x0_new,y0_new,file=files)
    print('saving new x0y0 message to x0y0_new.txt')
with open('parameters_new.txt','w+') as files:
    E0,e0,Waz,Wze=asa.run1(ra1,dec1,x,y,x0,y0,time0,loc,f)
    for i in range(20):
        E00=E0
        E0,a0,e0=asa.iteration(x,y,x0_new,y0_new,Waz,Wze,E0,e0,'iter')
        #print(i+1,'E=%.05f, a0=%.05f, e=%.05f, dE:%.03f%%'%(E0,a0,e0,(E00-E0)/E00*100))
    Caz,Cze,E,a0,e,k1,k2,k3,k4=asa.iteration(x,y,x0_new,y0_new,Waz,Wze,E0,e0,'')
    print(E,a0,e,k1,k2,k3,k4,file=files)
    print('saving new parameters to parameter_new.txt')


