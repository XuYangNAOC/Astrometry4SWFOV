from astropy.io import fits
import numpy as np
import math
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from scipy.optimize import curve_fit
import re
import all_sky_astro as asa

### step 2, cross match catalog
## get basic info
# parameters
x0,y0=np.loadtxt('x0y0.txt',dtype=float)
E,a0,e,k1,k2,k3,k4=np.loadtxt('parameters.txt',dtype=float)

# time
time0='2017-05-24 19:00:55'
time = Time(time0)

# location
latitude=38.330278
longitude=74.896667
height=4565
loc = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height*u.m)

# ccd info (for cross match ref)
ref_fwhm = 3            # the max_seperation = ref_fwhm*pix_scale*scale_unit
pix_scale = 8          #
scale_unit = 1*u.arcmin #

# target lists
reflist='tycho_v6_fullsky.txt'
targ0='XJCAM_20170524190039_1.cat'

# get time from file name (or just give a parameter as above)
#targ1=targ0.split('_')
#targ2=re.split('(..)',targ0.split('_')[1])
#time0=targ2[1]+targ2[3]+'-'+targ2[5]+'-'+targ2[7]+' '+targ2[9]+':'+targ2[11]+':'+str(int(targ2[13])+17)

## read catalog
# read target catalog
hdul = fits.open(targ0)

# check basic info of fits file
#hdul.info()        
#col=hdul[2].columns
#print(col)

# read data
data=hdul[2].data
x=data['XWIN_IMAGE']
y=data['YWIN_IMAGE']
mag=data['MAG_APER'].T[1]   # for multi apertures
#mag=data['MAG_APER']       # for single aperture
fwhm=data['FWHM_IMAGE']

# filter
r=np.sqrt((x-x0)**2+(y-y0)**2)
index1=set(list(np.where(mag>5.5)[0]))
index2=set(list(np.where(fwhm<0.1)[0]))
index3=set(list(np.where(fwhm>10)[0]))
index4=set(list(np.where(r>700)[0]))
#index5=set(list(np.where(mag<1)[0]))
index=np.array(list(index1|index2|index3|index4))#|index5))

x=np.delete(x,index)
y=np.delete(y,index)
r=np.delete(r,index)
mag=np.delete(mag,index)

# get ra dec from calculation
ra,dec,Caz,Cze,Cu,Cb=asa.xytoradec(x,y,x0,y0,E,a0,e,k1,k2,k3,k4,time0,loc)

## cross match
# read reference catalog
Rra,Rdec,Rmag=np.loadtxt(reflist,skiprows=1,dtype=float,unpack=True)

# get wcs
Rwcs=SkyCoord(ra=Rra*u.degree,dec=Rdec*u.degree)    # R for Reference
Cwcs=SkyCoord(ra=ra*u.degree,dec=dec*u.degree)      # C for Calculation
idx,d2d,d3d=Cwcs.match_to_catalog_3d(Rwcs)

max_sep=ref_fwhm*pix_scale*u.arcmin
sep_c=d2d<max_sep
Cmatch=Cwcs[sep_c]
Rmatch=Rwcs[idx[sep_c]]

cm=mag[sep_c]
rm=Rmag[idx[sep_c]]

cx=x[sep_c]
cy=y[sep_c]
cr=r[sep_c]

# save result
f=open('ref_'+targ0+'.txt','w+')
for i in range(len(cx)):
    print(i,Rmatch.ra.degree[i],Rmatch.dec.degree[i],cx[i],cy[i],file=f)
f.close()

## plot chec
"""
import matplotlib.pyplot as plt
# check Ref projection
Raz,Rze,Ru,Rb=asa.radectoub(Rra,Rdec,x0,y0,E,a0,e,k1,k2,k3,k4,time0)
cb_az=Cb[sep_c]
cu_ze=Cu[sep_c]
rb_az=Rb[idx[sep_c]]
ru_ze=Ru[idx[sep_c]]

c1=cb_az
c2=cu_ze
w1=rb_az
w2=ru_ze
"""
