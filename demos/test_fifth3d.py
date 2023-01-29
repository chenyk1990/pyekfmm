## This DEMO is a 3D example [x,y,z] with varying velocity and with multi-shots 
#  This example is benchmarked with pure C-code version to be correct.
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import time


z=np.linspace(0,2,201);
v=1.5+0.2*z;
vel=np.zeros([201,201,201],dtype='float32');
for i2 in range(201):
	for i3 in range(201):
		vel[:,i2,i3]=v;
vel=np.transpose(vel,[2,1,0]); ##from zxy to xyz 
vel=vel.reshape(201*201*201,1,order='F');
#see print(vel[0,0,:])

shot1=np.array([0.0,0,0],dtype='float32');
shot2=np.array([0.5,0,0],dtype='float32');
shot3=np.array([1.0,0,0],dtype='float32');
shots=np.concatenate([shot1,shot2,shot3],axis=0).reshape(3,3)

tic = time.perf_counter()
t=fmm.eikonal(vel,xyz=shot1,ax=[0,0.01,201],ay=[0,0.01,201],az=[0,0.01,201],order=2);
toc = time.perf_counter()
print(f"C version takes {toc - tic:0.4f} seconds");
time=t.reshape(201,201,201,1,order='F'); #[x,y,z]

import matplotlib.pyplot as plt
plt.imshow(time[:,50,:,0].transpose()); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Z');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()


## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',2.03724,0,0.329008,0.108246])
