import pyekfmm as fmm
import numpy as np
# t=fmm.eikonal(1,1);

vel=3.0*np.ones([501*501,1],dtype='float32');
t=fmm.eikonal(vel);




time=t.reshape(501,501,order='F');

import matplotlib.pyplot as plt

plt.imshow(time);
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()

