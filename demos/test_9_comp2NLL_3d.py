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

## axes information
nx=181;
ny=181;
nz=25;

dx=1;
dy=1;
dz=1;

ox=-90
oy=-90
oz=-3;

# number of events (calculations)
ne=1;

## lonlat corresponding to [x=0,y=0] 
LatOrig=39.665000;
LongOrig=-119.690000; #lon is x; lat is y

## station
stx=-57.291076;
sty=22.598270;
stz=-1.743000;

## reference time from NLL
file_time='/Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-grd3/g2t3.P.BEK.time.buf'
fd = open(file_time,'rb')
data = np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny],order='F')    #[zxy]

## 3D vel from NLL
file_vel='/Users/chenyk/softs/GrowClust3D.jl/examples/data/nll/nll-run/v2g3.inp.P.mod.buf'
fd = open(file_vel,'rb')
vel = 1./np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny],order='F')  #[zxy]

## plot velocity
plt.figure();
plt.imshow(vel[:,:,10],aspect='auto');
plt.jet();
plt.colorbar()
plt.title('Velocity model')
plt.show();


## pyekfmm
vel=np.swapaxes(vel,0,2).reshape([nz*nx*ny,1],order='F'); #transpose to [xyz]

shot=np.array([stx,sty,stz],dtype='float32');
t=fmm.eikonal(vel,xyz=shot,ax=[ox,dx,nx],ay=[oy,dy,ny],az=[oz,dz,nz],order=1);
time=t.reshape(nx,ny,nz,order='F'); #[x,y,z]
time=np.swapaxes(time,0,2); #[z,x,y]


## Plot
plt.figure();
plt.imshow(data[:,:,10],aspect='auto');
plt.jet();
plt.title('NLL')
plt.show();

plt.figure();
plt.imshow(time[:,:,10],aspect='auto');
plt.jet();
plt.title('pyekfmm')
plt.show();

import matplotlib.pyplot as plt
plt.figure();
plt.imshow(np.concatenate((data[:,:,10],time[:,:,10],data[:,:,10]-time[:,:,10]),axis=1),aspect='auto');
plt.jet();
plt.title('NLL V.S. pyekfmm V.S. Difference')
plt.colorbar();
plt.show();

dif=np.abs(data-time);
print("EKFMM and NLL difference:")
print("First order: std: %g, max:%g, min:%g, mean:%g, norm:%g \n"%(dif.std(),dif.max(),dif.min(),dif.mean(),np.linalg.norm(dif)))


## second order
t=fmm.eikonal(vel,xyz=shot,ax=[ox,dx,nx],ay=[oy,dy,ny],az=[oz,dz,nz],order=2);
time=t.reshape(nx,ny,nz,order='F'); #[x,y,z]
time=np.swapaxes(time,0,2); #[z,x,y]

dif=np.abs(data-time);
print("Second order: std: %g, max:%g, min:%g, mean:%g, norm:%g \n"%(dif.std(),dif.max(),dif.min(),dif.mean(),np.linalg.norm(dif)))


## Conclusion:
# first order difference is closer, meaning that NLL is up to the first-order accuracy
# while pyekfmm can have second-order accuracy (more accurate)
# 










