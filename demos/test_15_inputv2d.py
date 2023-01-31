## This DEMO is a 2D example [x,z] with input velocity and with one shot and one station and output the ray path
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt

## download this model from 
# https://github.com/aaspip/data/blob/main/vsmooth5063.bin

## input parameter (please note the axis sequence, input: zx or zxy)
file_time='./vsmooth5063.bin'
fd = open(file_time,'rb')
vel = np.fromfile(fd, dtype = np.float32).reshape([50,63],order='F')    #[zxy]
vel1d=vel.transpose().flatten(order='F');#first axis (vertical) is x, second is z, required by pyekfmm


## specify parameters
deg2km=6371*2*np.pi/360.0
x0=109*deg2km;
dx=0.75*deg2km;
z0=-46*deg2km;
dz=0.75*deg2km;
nx=63;
nz=50;

## maximum lon and lat in km
xm=x0+dx*nx;
zm=z0+dz*nz;

## Event location in km
xe=x0+dx*10; #lon
ze=z0+dz*10; #lat

## Event location in lat/lon
elon=146.66660;
elat=-42.30800;
## Station location in lat/lon
slon=147.32040;#This one is too close (you can uncomment it to see how it goes)
slat=-42.90990;#This one is too close (you can uncomment it to see how it goes)

slon=115.92670;#This one is further
slat=-31.98350;#This one is further
   
slon=128.92670;#A fake one
slat=-18.98350;#A fake one

## Event and Station in km
xe=elon*deg2km;
ze=elat*deg2km;
xs=slon*deg2km;
zs=slat*deg2km;

## Traveltime calculation
t=fmm.eikonal(vel1d,xyz=np.array([xe,0,ze]),ax=[x0,dx,nx],ay=[0,0.01,1],az=[z0,dz,nz],order=2);
time=t.reshape(nx,nz,order='F');#first axis (vertical) is x, second is z
time=time.transpose(); #z,x

## Ray tracing begins
tz,tx = np.gradient(time)

## trace one ray (source and receiver locations are in decimal grids)
stationx=(slon-x0/deg2km)/0.75+1
stationz=(slat-z0/deg2km)/0.75+1
eventx=(elon-x0/deg2km)/0.75+1
eventz=(elat-z0/deg2km)/0.75+1

## key function: stream2d
paths,nrays=fmm.stream2d(-tx,-tz, stationx, stationz, step=0.1, maxvert=10000)
print('Before trim',paths.shape)
## trim the rays (remove redundant rays around source) and add the source point
paths=fmm.trimrays(paths,start_points=np.array([eventx,eventz]),T=0.2)
print('After trim',paths.shape)

## Visualization
fig = plt.figure(figsize=(16, 8))
ax=plt.subplot(1,2,1)
plt.imshow(vel,cmap=plt.cm.jet, interpolation='none', extent=[x0/deg2km,xm/deg2km,zm/deg2km,z0/deg2km]); 
plt.plot(xe/deg2km,ze/deg2km,'*r',markersize=10)
plt.plot(xs/deg2km,zs/deg2km,'vb',markersize=10)
plt.xlabel('Lon (deg)');plt.ylabel('Lat (deg)');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6,label='Velocity (km/s)');
# plot the ray
plt.plot((paths[0,:]-1)*dx/deg2km+x0/deg2km,(paths[1,:]-1)*dz/deg2km+z0/deg2km,'g--',markersize=20);
plt.gca().invert_yaxis()

ax=plt.subplot(1,2,2)
plt.imshow(time,cmap=plt.cm.jet, interpolation='none', extent=[x0/deg2km,xm/deg2km,zm/deg2km,z0/deg2km]);
plt.plot(xe/deg2km,ze/deg2km,'*r',markersize=10)
plt.plot(xs/deg2km,zs/deg2km,'vb',markersize=10)
plt.xlabel('Lon (deg)');plt.ylabel('Lat (deg)');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
# plot the ray
plt.plot((paths[0,:]-1)*dx/deg2km+x0/deg2km,(paths[1,:]-1)*dz/deg2km+z0/deg2km,'g--',markersize=20);
plt.gca().invert_yaxis()
plt.savefig('test_15_inputv2d.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.show()



