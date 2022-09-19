import pyekfmm as fmm
import numpy as np
# t=fmm.eikonal(1,1);

vel=3.0*np.ones([101*101*101,1],dtype='float32');

# vel=3.0*np.ones([501*501,1],dtype='float32');
# t=fmm.eikonal(vel,xyz=np.array([0,0,0]),ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2);

t=fmm.fmm.eikonal(vel,xyz=np.array([0,0,0]));


# time=t.reshape(501,501,1,order='F');
time=t.reshape(101,101,101,order='F');

import matplotlib.pyplot as plt

plt.imshow(time[:,:,50]);
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

