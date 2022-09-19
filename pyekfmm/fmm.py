def eikonal(vel,xyz,ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2):
	'''
	EIKONAL: Fast marching eikonal solver (3-D)
	
	INPUT
	vel: 1D numpy array (nx*ny*nz)
	xyz: 1D/2D numpy array (one event: 1x3 or multi-event: ne x 3)
	ax: axis x [ox,dx,nx]
	ay: axis y [oy,dy,ny]
	az: axis z [oz,dz,nz]
	
	OUTPUT
	times: traveltime in xyz respectively (1D numpy array)
		   (one event: nx*ny*nz or multi-event: nx*ny*nz*ne)
	
	EXAMPLE
	demos/test_first.py
	
	COPYRIGHT
	Yangkang Chen, 2022, The University of Texas at Austin
	
	MODIFICATIONS
	[1] By Yangkang Chen, Sep, 2022
	
	'''
	import numpy as np
	if xyz.size == 3:
		from eikonalc import eikonalc_oneshot
		x=xyz[0];y=xyz[1];z=xyz[2];
		times=eikonalc_oneshot(vel,x,y,z,ax[0],ay[0],az[0],ax[1],ay[1],az[1],ax[2],ay[2],az[2],order);
	else:
		from eikonalc import eikonalc_multishots
		[ne,ndim]=xyz.shape;#ndim must be 3
		x=xyz[:,0];y=xyz[:,1];z=xyz[:,2];
		x=np.expand_dims(x,1);
		y=np.expand_dims(y,1);
		z=np.expand_dims(z,1);
		times=eikonalc_multishots(vel,x,y,z,ax[0],ay[0],az[0],ax[1],ay[1],az[1],ax[2],ay[2],az[2],order);
		
	return times
	
	
def eikonal_surf(vel,xyz,ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2):
	'''
	EIKONAL_SURF: Fast marching eikonal solver (3-D) and recording the traveltimes on the surface.
	
	INPUT
	vel: 1D numpy array (nx*ny*nz)
	xyz: 1D/2D numpy array (one event: 1x3 or multi-event: ne x 3)
	ax: axis x [ox,dx,nx]
	ay: axis y [oy,dy,ny]
	az: axis z [oz,dz,nz]
	
	OUTPUT
	times: traveltime in xy of nshots (1D numpy array)
		   size of times is nx*ny*nshots
	
	EXAMPLE
	demos/test_xxx.py
	
	COPYRIGHT
	Yangkang Chen, 2022, The University of Texas at Austin
	
	MODIFICATIONS
	[1] By Yangkang Chen, Sep, 2022
	'''
	from eikonalc import eikonalc_surf
	[ne,ndim]=xyz.shape;#ndim must be 3
	x=xyz[:,0];y=xyz[:,1];z=xyz[:,2];
	times=eikonalc_surf(vel,x,y,z,ax[0],ay[0],az[0],ax[1],ay[1],az[1],ax[2],ay[2],az[2],order);
		
	return times
	
