## This DEMO is a 2D example [t,p] (lat,lon) for tomo test
# 
#  COPYRIGHT: Yunfeng Chen, 2023

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from obspy.geodetics.base import locations2degrees   


## download this model from 
# https://github.com/aaspip/data/blob/main/vsmooth5063.bin

## input parameter (please note the axis sequence, input: zx or zxy)
# file_time='./vsmooth2131.bin'
# fd = open(file_time,'rb')
# vel = np.fromfile(fd, dtype = np.float32).reshape([21,31],order='F')    #[zxy]
ddx=0.025;
ddz=0.025;
nx=int((0.25/ddx)*30+1);
nz=int((0.25/ddz)*20+1);
vel = np.ones((nz,nx),dtype = np.float32)*3. #[rtp]
vel1d=vel.flatten(order='F');
# [r,t,p], r is radius, t is lat, p is lon 
# in 2D surface case, first axis is lat, second axis is lon, required by pyekfmm

## specify parameters
deg2km=6371*2*np.pi/360.0
x0=119;
dx=ddx;
z0=-35;
dz=ddz;
# nx=31;
# nz=21;

x=[(x0+i*dx) for i in range(nx)]
z=[(z0+i*dz) for i in range(nz)]

## maximum lon and lat in km
xm=x0+dx*nx;
zm=z0+dz*nz;

## Read in event and station locations
pairs = np.loadtxt('travel_time_AL.dat')
fid=open('raypath_AL.dat','w')
fid1=open('travel_time_syn_AL.dat','w')
distance=[]
travel_time=[]
# for i,p in enumerate(pairs[:1,:]):
for i,p in enumerate(pairs):

    if i % 500 == 0:
        print(f'tracing the {i} ray')
    ## Event location in lat/lon
    elon=p[0];
    elat=p[1];
    ## Station location in lat/lon  
    slon=p[2]
    slat=p[3]

    if elat>-30.5 or elat<-34 or elon>126 or elon<120 or \
        slat>-30.5 or slat<-34 or slon>126 or slon<120:
        continue
    ## Event and Station in km
    xe=elon;
    ze=elat;
    xs=slon;
    zs=slat;
    
    distance.append(locations2degrees(elat,elon,slat,slon)*deg2km)
    ## Traveltime calculation (I preserve the cartesian version for comparison)
#     t=fmm.eikonal(vel1d,xyz=np.array([xe,0,ze]),ax=[x0,dx,nx],ay=[0,0.01,1],az=[z0,dz,nz],order=2);
    t=fmm.eikonal_rtp(vel1d,rtp=np.array([6371,ze+90,xe+180]),ar=[6371,6371,1],at=[z0+90,dz,nz],ap=[x0+180,dx,nx],order=2); 
    #Remember, in the spherical coordinates, t(lat):[0,180], p(lon):[0:360]
    time=t.reshape(nz,nx,order='F');#[r,t,p], r is radius, t is lat, p is lon
    
    # Extract the travel time at the station location
    f=RectBivariateSpline(z,x,time)
    tt=f(zs,xs).flatten()[0]
    
    travel_time.append(tt)
    ## Ray tracing begins
    tz,tx = np.gradient(time)

    ## trace one ray (source and receiver locations are in decimal grids)
    stationx=(slon-x0)/ddx+1
    stationz=(slat-z0)/ddz+1
    eventx=(elon-x0)/ddx+1
    eventz=(elat-z0)/ddz+1

    ## key function: stream2d
    paths,nrays=fmm.stream2d(-tx,-tz, stationx, stationz, step=0.1, maxvert=10000)
    ## trim the rays (remove redundant rays around source) and add the source point
    paths=fmm.trimrays(paths,start_points=np.array([eventx,eventz]),T=0.2)
    
    # Save the ray paths and travel times as text files
    for px,pz in zip(paths[0],paths[1]):
        px=(px-1)*dx+x0
        pz=(pz-1)*dz+z0
        fid.write('{0:f} {1:f}\n'.format(px,pz))
    fid.write('{0:f} {1:f}\n'.format(np.nan,np.nan))   
    fid1.write('{0:f} {1:f} {2:f} {3:f} {4:f}\n'.format(elon,elat,slon,slat,tt))
    
    ## Visualization
    if i == 0:
        fig = plt.figure(figsize=(6, 6))
        ax=plt.subplot(1,1,1)
        plt.imshow(vel,cmap=plt.cm.jet, interpolation='none', extent=[x0,xm,zm,z0]); 
        ax.set_xlabel('Lon (deg)');
        ax.set_ylabel('Lat (deg)');
        plt.colorbar(orientation='horizontal',shrink=0.6,label='Velocity (km/s)');
        plt.gca().invert_yaxis()
    # mark event and station locations
    # ax1.plot(xe,ze,'*r',markersize=10)
    # ax1.plot(xs,zs,'vb',markersize=10)
    # plot the ray
    ax.plot((paths[0,:]-1)*dx+x0,(paths[1,:]-1)*dz+z0,
              'k',linewidth=0.1,alpha=0.4);
# plt.savefig('test_17.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.show()
fid.close()
fid1.close()

# Plot distance-travel time curve
plt.figure()
plt.plot(distance, travel_time,'.')
plt.xlabel('Distance (km)')
plt.ylabel('Travel time (s)')
plt.savefig('curve2_sph.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.show()

