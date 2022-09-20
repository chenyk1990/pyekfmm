## This DEMO is a 3D example [x,y,z] with constant velocity and with multi-shots with surface recording only
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');

shot1=np.array([0.0,0,0],dtype='float32');
shot2=np.array([0.5,0,0],dtype='float32');
shot3=np.array([1.0,0,0],dtype='float32');
shots=np.concatenate([shot1,shot2,shot3],axis=0).reshape(3,3)

t=fmm.eikonal_surf(vel,xyz=shots,ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2);
time=t.reshape(101,101,3,order='F'); #[x,y]


import matplotlib.pyplot as plt
plt.imshow(time[:,:,0].transpose()); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Y');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

import matplotlib.pyplot as plt
plt.imshow(time[:,:,1].transpose()); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Y');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

import matplotlib.pyplot as plt
plt.imshow(time[:,:,2].transpose()); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Y');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',0.47041193, 0.0, 0.096254215, 0.009264874])