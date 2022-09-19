# import pyekfmm as fmm
from eikonalc import eikonalc

def eikonal(vel,xyz,ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2):
	'''
	eikonal: Fast marching eikonal solver (3-D)
	
	INPUT
	xyz: np.ndarray
	ax: axis x
	ay: axis y
	az: axis z
	
	OUTPUT
	times: traveltime (numpy array)
	
	EXAMPLE
	
	'''
# 	[ne,ndim]=xyz.shape;
	
	x=xyz[0];y=xyz[1];z=xyz[2];
	print('x=',x);
# 	times=eikonalc(vel,x,y,z,ax[0],ay[0],az[0],ax[1],ay[1],az[1],ax[2],ay[2],az[2],order);
	times=eikonalc(vel,x,y,z,ax[0],ay[0],az[0],ax[1],ay[1],az[1],ax[2],ay[2],az[2],order);
# 	times=eikonalc(vel);
# 	times=[]

	return times