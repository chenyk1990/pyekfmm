## This DEMO is a 2D example [x,z] with input velocity and with one shot and one station and output the ray path
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

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

x=[(x0+i*dx) for i in range(nx)]
z=[(z0+i*dz) for i in range(nz)]


## maximum lon and lat in km
xm=x0+dx*nx;
zm=z0+dz*nz;

## Read in event and station locations
pairs = np.loadtxt('travel_time.dat')
fid=open('raypath.dat','w')
fid1=open('travel_time_syn.dat','w')
for i,p in enumerate(pairs):
    if i % 500 == 0:
        print(f'tracing the {i} ray')
    ## Event location in lat/lon
    elon=p[0];
    elat=p[1];
    ## Station location in lat/lon  
    slon=p[2]
    slat=p[3]

    ## Event and Station in km
    xe=elon*deg2km;
    ze=elat*deg2km;
    xs=slon*deg2km;
    zs=slat*deg2km;

    ## Traveltime calculation
    t=fmm.eikonal(vel1d,xyz=np.array([xe,0,ze]),ax=[x0,dx,nx],ay=[0,0.01,1],az=[z0,dz,nz],order=2);
    time=t.reshape(nx,nz,order='F');#first axis (vertical) is x, second is z
    time=time.transpose(); #z,x
    
    # Extract the travel time at the station location
    f=RectBivariateSpline(z,x,time)
    tt=f(zs,xs).flatten()[0]
    tt1=p[4]
    
    ## Ray tracing begins
    tz,tx = np.gradient(time)

    ## trace one ray (source and receiver locations are in decimal grids)
    stationx=(slon-x0/deg2km)/0.75+1
    stationz=(slat-z0/deg2km)/0.75+1
    eventx=(elon-x0/deg2km)/0.75+1
    eventz=(elat-z0/deg2km)/0.75+1

    ## key function: stream2d
    paths,nrays=fmm.stream2d(-tx,-tz, stationx, stationz, step=0.1, maxvert=10000)
    ## trim the rays (remove redundant rays around source) and add the source point
    paths=fmm.trimrays(paths,start_points=np.array([eventx,eventz]),T=0.2)
    
    # Save the ray paths and travel times as text files
    for px,pz in zip(paths[0],paths[1]):
        px=(px-1)*dx/deg2km+x0/deg2km
        pz=(pz-1)*dz/deg2km+z0/deg2km
        fid.write('{0:f} {1:f}\n'.format(px,pz))
    fid.write('{0:f} {1:f}\n'.format(np.nan,np.nan))   
    fid1.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f}\n'.format(elon,elat,slon,slat,tt,tt1))
    
    ## Visualization
    if i == 0:
        fig = plt.figure(figsize=(6, 6))
        ax=plt.subplot(1,1,1)
        plt.imshow(vel,cmap=plt.cm.jet, interpolation='none', extent=[x0/deg2km,xm/deg2km,zm/deg2km,z0/deg2km]); 
        ax.set_xlabel('Lon (deg)');
        ax.set_ylabel('Lat (deg)');
        plt.colorbar(orientation='horizontal',shrink=0.6,label='Velocity (km/s)');
        plt.gca().invert_yaxis()
    # mark event and station locations
    # ax1.plot(xe/deg2km,ze/deg2km,'*r',markersize=10)
    # ax1.plot(xs/deg2km,zs/deg2km,'vb',markersize=10)
    # plot the ray
    ax.plot((paths[0,:]-1)*dx/deg2km+x0/deg2km,(paths[1,:]-1)*dz/deg2km+z0/deg2km,
              'k',linewidth=0.1,alpha=0.4);
plt.savefig('test_16_tomo_rays.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.show()
fid.close()
fid1.close()

