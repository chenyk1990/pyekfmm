
## This DEMO is a 3D example [x,y,z] with constant velocity and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');
t=fmm.eikonal(vel,xyz=np.array([0.5,0,0]),ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2);
time=t.reshape(101,101,101,order='F'); #[x,y,z]


import matplotlib.pyplot as plt
plt.imshow(time[:,50,:].transpose()); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Z');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()


