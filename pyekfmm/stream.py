import numpy as np

def stream2d(u,v, sx, sy, step=0.1, maxvert=10000):
	"""
	stream2d: draw stream lines along the steepest descent direction
	
	INPUT
	u:   	derivative of traveltime in x
	v:   	derivative of traveltime in z
	
	OUTPUT  
	
	Copyright (C) 2023 The University of Texas at Austin
	Copyright (C) 2023 Yangkang Chen
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	References   
	[1] Chen et al.
		
	 DEMO
	 demos/test_xxx.py
	"""
	
	xSize=u.shape[1]
	ySize=u.shape[0]
	
	[verts, numverts] = traceStreamUV  (u.flatten(order='F'),v.flatten(order='F'),xSize, ySize, sx, sy, step, maxvert);
	verts=verts.reshape(2,numverts,order='F');
	
	return verts,numverts
	
def traceStreamUV (ugrid, vgrid, xdim, ydim, sx, sy, step, maxvert):
	"""
	traceStreamUV: key programs
	
	INPUT
	ugrid:   	derivative of traveltime in x
	vgrid:   	derivative of traveltime in z
	
	OUTPUT  
	
	Copyright (C) 2023 The University of Texas at Austin
	Copyright (C) 2023 Yangkang Chen
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details: http://www.gnu.org/licenses/
	
	References   
	[1] Chen et al.
		
	 DEMO
	 demos/test_xxx.py
	"""
	numverts=0;
	x=sx-1;y=sy-1;
	verts=np.zeros([2*maxvert,1])
	while 1:
		if (x<0 or x>xdim-1 or y<0 or y>ydim-1 or numverts>=maxvert) : 
# 			print("First break");
			break;
		
		ix=int(np.floor(x))
		iy=int(np.floor(y))
		
		if ix == xdim-1:
			ix=ix-1;
			
		if iy == ydim-1:
			iy=iy-1;
			
		xfrac=x-ix;
		yfrac=y-iy;
		
		#weights for linear interpolation
		a=(1-xfrac)*(1-yfrac);
		b=(  xfrac)*(1-yfrac);
		c=(1-xfrac)*(  yfrac);
		d=(  xfrac)*(  yfrac);
		
		verts[2*numverts + 0] = x+1;
		verts[2*numverts + 1] = y+1;
		
		#if already been here, done
		if numverts>=2:
			if verts[2*numverts] == verts[2*(numverts-2)] and verts[2*numverts+1] == verts[2*(numverts-2)+1]:
				numverts=numverts+1;
# 				print("Second break");
				break;
		
		numverts=numverts+1;
		ui = ugrid[iy  +ydim*ix]*a + ugrid[iy  +ydim*(ix+1)]*b + ugrid[iy+1+ydim*ix]*c + ugrid[iy+1+ydim*(ix+1)]*d;
		vi = vgrid[iy  +ydim*ix]*a + vgrid[iy  +ydim*(ix+1)]*b + vgrid[iy+1+ydim*ix]*c + vgrid[iy+1+ydim*(ix+1)]*d;
		
		#calculate step size, if 0, done
		if abs(ui) > abs(vi):
			imax=abs(ui);
		else:
			imax=abs(vi);
			
		if imax==0:
			print("Third break");
			break;
		
		imax=step/imax;
		ui = ui*imax;
		vi = vi*imax;
		
		#update the current position
		x = x+ui;
		y = y+vi;
	
# 	print('numverts',numverts)
	verts=verts[0:2*numverts]
	
	return verts,numverts
	
def trimrays(paths, start_points, T=None):
	"""
	trimrays: trim rays (remove very close ray points around the source)
	
	%paths: [2 x ngrid]
	%start_points: [2 x 1], e.g., start_points=np.array([1,1])
	%T: threshold, e.g., dx,dz
	
	"""
	ngrid=paths.shape[1]
	start_points=np.expand_dims(start_points,1);
	
	d=np.sqrt(np.sum(np.power(paths-np.repeat(start_points,ngrid,axis=1),2),0));
	
	if T is None:
		T = d.max()/300/300;
		
	I=np.where(d<T)
	
	paths=np.delete(paths,I,1)
	
	paths=np.concatenate((paths,start_points),axis=1)
	
	return paths
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	