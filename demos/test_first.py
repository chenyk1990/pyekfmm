## This DEMO is a 2D example [x,z] with constant velocity and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([501*501,1],dtype='float32'); #velocity axis must be x,y,z respectively
t=fmm.eikonal(vel,xyz=np.array([2.5,0,0]),ax=[0,0.01,501],ay=[0,0.01,1],az=[0,0.01,501],order=2);
time=t.reshape(501,501,order='F');#first axis (vertical) is x, second is z


import matplotlib.pyplot as plt
plt.imshow(time.transpose()); #transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Z');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',1.8627454, 0.0, 0.42475563, 0.18041734])

