## This DEMO is a 2D example [x,z] with constant velocity in VTI media and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

# vel=3.80395*np.ones([501*501,1],dtype='float32'); #velocity axis must be x,y,z respectively
# velx=3.09354*np.ones([501*501,1],dtype='float32'); #velocity axis must be x,y,z respectively
# eta=0.340859*np.ones([501*501,1],dtype='float32'); #velocity axis must be x,y,z respectively
velz=3.09354*np.ones([201*201,1],dtype='float32'); #velocity axis must be x,y,z respectively
velx=3.80395*np.ones([201*201,1],dtype='float32'); #velocity axis must be x,y,z respectively
eta=0.340859*np.ones([201*201,1],dtype='float32'); #velocity axis must be x,y,z respectively
t=fmm.eikonalvti(velx,velz,eta,xyz=np.array([5,0,5]),ax=[0,0.05,201],ay=[0,0.05,1],az=[0,0.05,201],order=2);
time=t.reshape(201,201,order='F');#first axis (vertical) is x, second is z


import matplotlib.pyplot as plt
plt.imshow(time.transpose(),cmap=plt.cm.jet, interpolation='none', extent=[0,10,10,0]); #transpose so that first axis is z, second is x
plt.plot(5,5,'*r',markersize=10)
plt.xlabel('X (km)');plt.ylabel('Z (km)');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
plt.savefig('test_12_constvti.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.show()

## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',2.1514063, 0.0, 0.44167396, 0.19507588])

