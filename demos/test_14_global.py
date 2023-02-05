## This DEMO is a 3D spherical example [r,t,p] with constant velocity and with one shot
#  in a global scale 
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

##Event lat/lon/dep
evlat=-25+90
evlon=0+180
evdep=0

## size of rtp axes (number of samples)
numr=1
numt=1000
nump=2000

## beg/end in r
begr=6400
endr=6400

begt=0
endt=180

## beg/end in p
begp=0
endp=360

## sampling in r
dr=(endr-begr)/numr
dt=(endt-begt)/numt
dp=(endp-begp)/nump


vel=6.0*np.ones([numr*numt*nump,1],dtype='float32');
t=fmm.eikonal_rtp(vel,rtp=np.array([6400-evdep,evlat,evlon]),ar=[begr,dr,numr],at=[begt,dt,numt],ap=[begp,dp,nump],order=1);#spherical needs order=1
time=t.reshape(numr,numt,nump,order='F'); #[r,t,p]

import matplotlib.pyplot as plt
# plt.figure(figsize=(12, 7))
plt.imshow(time[numr-1,:,:],extent=[begp-180,endp-180,endt-90,begt-90],aspect='auto')#,clim=(0, 2000.0)); #rtp->tp
plt.gca().invert_yaxis()
plt.xlabel('Lon (deg)'); 
plt.ylabel('Lat (deg)');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
plt.plot(evlon-180,evlat-90,'*',color='r',markersize=15)
# plt.gca().invert_yaxis()
plt.savefig('test_14_global.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.savefig('test_14_global.pdf',format='pdf',dpi=300,bbox_inches='tight', pad_inches=0)
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',2140.937, 0.016649984, 75.647644, 5722.5664])
#print(['Correct result:',285.1792, 0.016649984, 75.70713, 5731.57]) #If using order=2



