## This DEMO is a 3D example [x,y,z] with varying velocity and with multi-shots 
#  This example is benchmarked with pure C-code version to be correct.
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import time


for ngrid in [2,4,8,16,32,64,128,256,512]:

	z=np.linspace(0,2,ngrid);
	v=1.5+0.2*z;
	vel=np.zeros([ngrid,ngrid,ngrid],dtype='float32');
	for i2 in range(ngrid):
		for i3 in range(ngrid):
			vel[:,i2,i3]=v;
	vel=np.transpose(vel,[2,1,0]); ##from zxy to xyz 
	vel=vel.reshape(ngrid*ngrid*ngrid,1,order='F');
#see print(vel[0,0,:])

	shot1=np.array([ngrid/2*0.01,ngrid/2*0.01,ngrid/2*0.01],dtype='float32');
	print('ngrid=',ngrid,'shot=',ngrid/2*0.01,ngrid/2*0.01,ngrid/2*0.01)
	tic = time.perf_counter()
	t=fmm.eikonal(vel,xyz=shot1,ax=[0,0.01,ngrid],ay=[0,0.01,ngrid],az=[0,0.01,ngrid],order=1);
	toc = time.perf_counter()
	print(f"{ngrid:2d} takes {toc - tic:0.8f} seconds");
# 	time=t.reshape(ngrid,ngrid,ngrid,1,order='F'); #[x,y,z]
# ngrid= 2 shot= 0.01 0.01 0.01
#  2 takes 0.00001112 seconds
# ngrid= 4 shot= 0.02 0.02 0.02
#  4 takes 0.00010137 seconds
# ngrid= 8 shot= 0.04 0.04 0.04
#  8 takes 0.00021634 seconds
# ngrid= 16 shot= 0.08 0.08 0.08
# 16 takes 0.00108170 seconds
# ngrid= 32 shot= 0.16 0.16 0.16
# 32 takes 0.00581083 seconds
# ngrid= 64 shot= 0.32 0.32 0.32
# 64 takes 0.06067825 seconds
# ngrid= 128 shot= 0.64 0.64 0.64
# 128 takes 0.75102152 seconds
# ngrid= 256 shot= 1.28 1.28 1.28
# 256 takes 10.77394834 seconds
# ngrid= 512 shot= 2.56 2.56 2.56
# 512 takes 147.03446876 seconds