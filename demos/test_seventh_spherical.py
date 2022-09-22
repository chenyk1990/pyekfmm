## This DEMO is a 3D spherical example [r,t,p] with constant velocity and with one shot
#  in a regional scale (e.g., Texas)
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

##Event lat/lon/dep
evlat=25
evlon=-108
evdep=100

## size of rtp axes (number of samples)
numr=1001
numt=200
nump=200

## beg/end in r
begr=6200
endr=6400

## beg/end in t
begt=20
endt=40

## beg/end in p
begp=-110
endp=-100

## sampling in r
dr=(endr-begr)/numr
dt=(endt-begt)/numt
dp=(endp-begp)/nump


vel=6.0*np.ones([numr*numt*nump,1],dtype='float32');
t=fmm.eikonal_rtp(vel,rtp=np.array([6400-evdep,evlat,evlon]),ar=[begr,dr,numr],at=[begt,dt,numt],ap=[begp,dp,nump],order=1);#spherical needs order=1
time=t.reshape(numr,numt,nump,order='F'); #[r,t,p]

import matplotlib.pyplot as plt
plt.imshow(time[numr-1,:,:],extent=[begp,endp,endt,begt],aspect='auto')#,clim=(0, 2000.0)); #rtp->tp
plt.gca().invert_yaxis()
plt.xlabel('Lon'); 
plt.ylabel('Lat');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.plot(evlon,evlat,'*',color='r',markersize=15)
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',286.06723, 0.016649984, 75.647644, 5722.5664])
#print(['Correct result:',285.1792, 0.016649984, 75.70713, 5731.57]) #If using order=2



