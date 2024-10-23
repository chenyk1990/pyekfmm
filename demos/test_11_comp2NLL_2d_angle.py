# Compare between PYEKFMM and NLL
# 
# Paths:
# 
# Station file
# /Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-run/nll-sta.txt 
#  GTSRCE   BEK LATLON 39.866600 -120.359600 0.000 1.743
#
# V file
# /Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-run/nll-vz.txt
# /Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-run/v2g3.inp.P.mod.buf
#
# Time file
# /Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-grd3/g2t3.P.BEK.time.buf
# /Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-grd3/g2t3.P.BEK.time.hdr

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt
import os

## axes information
nx=401;
ny=1;
nz=106;

dx=1;
dy=1;
dz=1;

ox=0
oy=0
oz=-5;

# number of events (calculations)
ne=1;

## lonlat corresponding to [x=0,y=0] 
LatOrig=61.00;
LongOrig=-150.00; #lon is x; lat is y

## station
stx=14.081350;
sty=9.933220;
stz=-0.390000;

## reference time from NLL
file_time=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.AK_RC01_--.time.buf'
fd = open(file_time,'rb')
time0=np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny*2],order='F')
time0=time0[:,:,0].reshape([nz,nx,ny],order='F')

file_time=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.AK_RC01_--.angle.buf'
fd = open(file_time,'rb')
data = np.fromfile(fd, dtype = np.uint16);#.reshape([nz,nx,ny*4],order='F')    #[zxy]
dip0=data[np.linspace(0,nz*nx*ny*2-2,nz*nx*ny,dtype='int')]/16/10;
azim0=data[np.linspace(1,nz*nx*ny*2-1,nz*nx*ny,dtype='int')]/10;

dip0=dip0.reshape([nz,nx,ny],order='F')
azim0=azim0.reshape([nz,nx,ny],order='F')

## 3D vel from NLL
file_vel=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/model/layer.P.mod.buf'
fd = open(file_vel,'rb')
vel = 1./np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny*2],order='F')  #[zxy]
vel=vel[:,:,0].reshape([nz,nx,ny],order='F')  #[zxy]

## plot time
# plt.figure();
# plt.imshow(time0[:,:,0],aspect='auto');
# plt.jet();
# plt.colorbar()
# plt.title('NLL Time')
# plt.show();
# 
# ## plot velocity
# plt.figure();
# plt.imshow(vel[:,:,0],aspect='auto');
# plt.jet();
# plt.colorbar()
# plt.title('Velocity model')
# plt.show();
# 
# ## plot dip
# plt.figure();
# plt.imshow(dip0[:,:,0],aspect='auto');
# plt.jet();
# plt.colorbar()
# plt.title('Dip model')
# plt.show();
# 
# ## plot azimuth
# plt.figure();
# plt.imshow(azim0[:,:,0],aspect='auto');
# plt.jet();
# plt.colorbar()
# plt.title('Azim model')
# plt.show();


## pyekfmm
vel=np.swapaxes(vel,0,2).reshape([nz*nx*ny,1],order='F'); #transpose to [xyz]
shot=np.array([0,0,0],dtype='float32');
# shot=np.array([stx,sty,stz],dtype='float32');
t,d,a=fmm.eikonal(vel,xyz=shot,ax=[ox,dx,nx],ay=[oy,dy,ny],az=[oz,dz,nz],order=1,angle=True);
time=t.reshape(nx,ny,nz,order='F'); #[x,y,z]
time=np.swapaxes(time,0,2); #[z,x,y]

dip=d.reshape(nx,ny,nz,order='F');
dip=np.swapaxes(dip,0,2); #[z,x,y]

azim=a.reshape(nx,ny,nz,order='F');
azim=np.swapaxes(dip,0,2); #[z,x,y]

print("Angle difference: ",np.linalg.norm(dip0[:,:,0]-dip[:,0,:]))
print("Time difference: ",np.linalg.norm(time0[:,:,0]-time[:,0,:]))

plt.figure();
plt.subplot(1,3,1)
plt.imshow(dip0[:,:,0],aspect='auto');
plt.jet();
plt.title('NLL')

plt.subplot(1,3,2)
plt.imshow(dip[:,0,:],aspect='auto');plt.gca().set_yticks([]);
plt.jet();
plt.title('pyekfmm')

plt.subplot(1,3,3)
plt.imshow(dip0[:,:,0]-dip[:,0,:],clim=[0,50],aspect='auto');plt.gca().set_yticks([]);
plt.jet();
plt.title('NLL-pyekfmm')
plt.colorbar(orientation='vertical',shrink=0.6);
plt.show();

plt.figure();
plt.subplot(1,3,1)
plt.imshow(azim0[:,:,0],clim=(0, 50),aspect='auto');
plt.jet();
plt.title('NLL')

plt.subplot(1,3,2)
plt.imshow(azim[:,0,:],clim=(0, 50),aspect='auto');plt.gca().set_yticks([]);
plt.jet();
plt.title('pyekfmm')

plt.subplot(1,3,3)
plt.imshow(dip0[:,:,0]-dip[:,0,:],clim=[0,50],aspect='auto');plt.gca().set_yticks([]);
plt.jet();
plt.title('NLL-pyekfmm')
plt.colorbar(orientation='vertical',shrink=0.6);
plt.show();


# 
# 
# ## Plot
# plt.figure();
# plt.imshow(data[:,:,10],aspect='auto');
# plt.jet();
# plt.title('NLL')
# plt.show();
# 
plt.figure();
plt.subplot(1,3,1)
plt.imshow(time0[:,:,0],aspect='auto');
plt.jet();
plt.title('NLL')

plt.subplot(1,3,2)
plt.imshow(time[:,0,:],aspect='auto');plt.gca().set_yticks([]);
plt.jet();
plt.title('pyekfmm')

plt.subplot(1,3,3)
plt.imshow(time0[:,:,0]-time[:,0,:],clim=[0,200],aspect='auto');plt.gca().set_yticks([]);
plt.jet();
plt.title('NLL-pyekfmm')
plt.colorbar(orientation='vertical',shrink=0.6);
plt.show();




# 
# import matplotlib.pyplot as plt
# plt.figure();
# plt.imshow(np.concatenate((data[:,:,10],time[:,:,10],data[:,:,10]-time[:,:,10]),axis=1),aspect='auto');
# plt.jet();
# plt.title('NLL V.S. pyekfmm V.S. Difference')
# plt.colorbar();
# plt.show();
# 
# dif=np.abs(data-time);
# print("EKFMM and NLL difference:")
# print("First order: std: %g, max:%g, min:%g, mean:%g, norm:%g \n"%(dif.std(),dif.max(),dif.min(),dif.mean(),np.linalg.norm(dif)))
# 
# 
# ## second order
# t=fmm.eikonal(vel,xyz=shot,ax=[ox,dx,nx],ay=[oy,dy,ny],az=[oz,dz,nz],order=2);
# time=t.reshape(nx,ny,nz,order='F'); #[x,y,z]
# time=np.swapaxes(time,0,2); #[z,x,y]
# 
# dif=np.abs(data-time);
# print("Second order: std: %g, max:%g, min:%g, mean:%g, norm:%g \n"%(dif.std(),dif.max(),dif.min(),dif.mean(),np.linalg.norm(dif)))


## Conclusion:
# first order difference is closer, meaning that NLL is up to the first-order accuracy
# while pyekfmm can have second-order accuracy (more accurate)
# 










