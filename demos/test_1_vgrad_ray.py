## This DEMO is a 2D example [x,z] with velocity gradient and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt

v1=1;
v2=3;
nz=501;
nx=501;
dx=0.01;
dz=0.01;
# vel=3.0*np.ones([501*501,1],dtype='float32'); #velocity axis must be x,y,z respectively
v=np.linspace(v1,v2,nz);
v=np.expand_dims(v,1);
h=np.ones([1,nx])
vel=np.multiply(v,h,dtype='float32'); #z,x
# plt.imshow(vel);plt.jet();plt.show()

t=fmm.eikonal(vel.transpose().flatten(order='F'),xyz=np.array([0,0,0]),ax=[0,dx,nx],ay=[0,0.01,1],az=[0,dz,nz],order=2);
time=t.reshape(nx,nz,order='F');#first axis (vertical) is x, second is z
time=time.transpose(); #z,x

# tz=np.gradient(time,axis=1);
# tx=np.gradient(time,axis=0);
# # or
tz,tx = np.gradient(time)




plt.figure();
plt.imshow(time,cmap=plt.cm.jet, interpolation='none', extent=[0,5,5,0]); #transpose so that first axis is z, second is x
plt.plot(0,0.0,'*r',markersize=20);
plt.xlabel('X (km)');plt.ylabel('Z (km)');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');

for ii in range(1,502,50):
	paths,nrays=fmm.stream2d(-tx,-tz, 501, ii, step=0.1, maxvert=10000)
	plt.plot(500*dx,(ii-1)*dz,'vb',markersize=15);
	## plot rays
	plt.plot((paths[0,:]-1)*dx,(paths[1,:]-1)*dz,'g--',markersize=20);

plt.savefig('test_1_vgrad_ray.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
plt.savefig('test_1_vgrad_ray.pdf',format='pdf',dpi=300,bbox_inches='tight', pad_inches=0)

plt.show()





## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',4.4039335, 0.0, 0.7773146, 0.604218])



