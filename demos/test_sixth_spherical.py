## This DEMO is a 3D spherical example [r,t,p] with constant velocity and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');
t=fmm.eikonal_rtp(vel,rtp=np.array([0.5,0,0]),ar=[0,0.01,101],at=[-90,1.8,101],ap=[-180,3.6,101],order=1);#spherical needs order=1
time=t.reshape(101,101,101,order='F'); #[r,t,p]

import matplotlib.pyplot as plt
plt.imshow(time[:,:,50])#,clim=(0, 1.0)); #rtp->rt
plt.xlabel('T');plt.ylabel('R');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',0.37473878, 1.9868216e-08, 0.06834167, 0.004670584])

## Below DEMO is a 3D spherical example [r,t,p] with constant velocity and with multi-shots
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');

shot1=np.array([0.0,0,0],dtype='float32');
shot2=np.array([0.5,0,0],dtype='float32');
shot3=np.array([1.0,0,0],dtype='float32');
shots=np.concatenate([shot1,shot2,shot3],axis=0).reshape(3,3)

t=fmm.eikonal_rtp(vel,rtp=shots,ar=[0,0.01,101],at=[-90,1.8,101],ap=[-180,3.6,101],order=1);#spherical needs order=1
time=t.reshape(101,101,101,3,order='F'); #[r,t,p]

import matplotlib.pyplot as plt
plt.imshow(time[:,:,50,0])#,clim=(0, 1.0)); #rtp->rt
plt.xlabel('T');plt.ylabel('R');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

import matplotlib.pyplot as plt
plt.imshow(time[:,:,50,1])#,clim=(0, 1.0)); #rtp->rt
plt.xlabel('T');plt.ylabel('R');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

import matplotlib.pyplot as plt
plt.imshow(time[:,:,50,2])#,clim=(0, 1.0)); #rtp->rt
plt.xlabel('T');plt.ylabel('R');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',0.4776695, 0.0, 0.09843297, 0.00968905])



