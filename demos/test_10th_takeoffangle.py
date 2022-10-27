## This DEMO is a 3D example [x,y,z] with constant velocity and with multi-shots 
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');

shot1=np.array([0.0,0,0],dtype='float32');
shot2=np.array([0.5,0,0],dtype='float32');
shot3=np.array([1.0,0,0],dtype='float32');
shots=np.concatenate([shot1,shot2,shot3],axis=0).reshape(3,3)

t,d,a=fmm.eikonal(vel,xyz=shots,ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2,angle=True);
time=t.reshape(101,101,101,3,order='F'); #[x,y,z]

## Takeoff angles
dip=d.reshape(101,101,101,3,order='F'); #[x,y,z]
azim=a.reshape(101,101,101,3,order='F'); #[x,y,z]

# import matplotlib.pyplot as plt
# plt.imshow(time[:,50,:,0].transpose()); #xyz->xz, transpose so that first axis is z, second is x
# plt.xlabel('X');plt.ylabel('Z');
# plt.jet()
# plt.colorbar(orientation='horizontal',shrink=0.6);
# plt.show()
# 
# import matplotlib.pyplot as plt
# plt.imshow(time[:,50,:,1].transpose()); #xyz->xz, transpose so that first axis is z, second is x
# plt.xlabel('X');plt.ylabel('Z');
# plt.jet()
# plt.colorbar(orientation='horizontal',shrink=0.6);
# plt.show()
# 
# import matplotlib.pyplot as plt
# plt.imshow(time[:,50,:,2].transpose()); #xyz->xz, transpose so that first axis is z, second is x
# plt.xlabel('X');plt.ylabel('Z');
# plt.jet()
# plt.colorbar(orientation='horizontal',shrink=0.6);
# plt.show()

## Take off angle
import matplotlib.pyplot as plt

## First
plt.subplot(1,3,1)
plt.imshow(time[:,50,:,0].transpose(),cmap='jet'); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Z');
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.subplot(1,3,2)
plt.imshow(dip[:,50,:,0].transpose(),cmap='jet',clim=(0, 150)); 
plt.xlabel('X');plt.ylabel('Z');plt.gca().set_yticks([]);
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.subplot(1,3,3)
plt.imshow(azim[:,50,:,0].transpose(),cmap='jet',clim=(0, 360)); 
plt.xlabel('X');plt.ylabel('Z');plt.gca().set_yticks([]);
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Second
plt.subplot(1,3,1)
plt.imshow(time[:,50,:,1].transpose(),cmap='jet'); 
plt.xlabel('X');plt.ylabel('Z');
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.subplot(1,3,2)
plt.imshow(dip[:,50,:,1].transpose(),cmap='jet',clim=(0, 150)); 
plt.xlabel('X');plt.ylabel('Z');plt.gca().set_yticks([]);
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.subplot(1,3,3)
plt.imshow(azim[:,50,:,1].transpose(),cmap='jet',clim=(0, 360)); 
plt.xlabel('X');plt.ylabel('Z');plt.gca().set_yticks([]);
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Third
plt.subplot(1,3,1)
plt.imshow(time[:,50,:,2].transpose(),cmap='jet'); 
plt.xlabel('X');plt.ylabel('Z');
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.subplot(1,3,2)
plt.imshow(dip[:,50,:,2].transpose(),cmap='jet',clim=(0, 150)); 
plt.xlabel('X');plt.ylabel('Z');plt.gca().set_yticks([]);
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.subplot(1,3,3)
plt.imshow(azim[:,50,:,2].transpose(),cmap='jet',clim=(0, 360)); 
plt.xlabel('X');plt.ylabel('Z');plt.gca().set_yticks([]);
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',0.5770362, 0.0, 0.094430424, 0.008917105])

