
import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt
import os

if os.path.isdir(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm') == False:  
	os.makedirs(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm',exist_ok=True)
	
## axes information
nx=401;
ny=1;
nz=106;

dx=1;
dy=1;
dz=1;

ox=0
oy=0
oz=-5;

# number of events (calculations)
ne=1;

## lonlat corresponding to [x=0,y=0] 
LatOrig=61.00;
LongOrig=-150.00; #lon is x; lat is y

## station
stx=14.081350;
sty=9.933220;
stz=-0.390000;


## 3D vel from NLL
file_vel=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/model/layer.P.mod.buf'
fd = open(file_vel,'rb')
vel = 1./np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny*2],order='F')  #[zxy]
vel=vel[:,:,0].reshape([nz,nx,ny],order='F')  #[zxy]
vel=np.swapaxes(vel,0,2).reshape([nz*nx*ny,1],order='F'); #transpose to [xyz]

## read station list
fid=open(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/obs/station_coordinates.txt');
lines=fid.readlines();
lines=lines[1:]
stnames=[ii.split()[1] for ii in lines];

os.system("cp -r /Users/chenyk/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.mod.* /Users/chenyk/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm/")
os.system("cp -r /Users/chenyk/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.*.angle.* /Users/chenyk/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm/")
for ii in range(len(stnames)):
	name=stnames[ii]
	os.system("cp -r /Users/chenyk/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.%s.time.hdr /Users/chenyk/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm/layer.P.%s.time.hdr"%(name,name))

	fid=open(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.%s.time.hdr'%name)
	line=fid.readlines()[1];
	shot=np.array(line.split()[1:4],dtype='float32')
	print(name,shot[0],shot[1],shot[2])
	
	fd=open(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.%s.time.buf'%name)
	time0=np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny*2],order='F')
	time0=time0[:,:,0].reshape([nz,nx,ny],order='F')
	
	file_time=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.%s.angle.buf'%name
	fd = open(file_time,'rb')
	data = np.fromfile(fd, dtype = np.uint16);#.reshape([nz,nx,ny*4],order='F')    #[zxy]
	dip0=data[np.linspace(0,nz*nx*ny*2-2,nz*nx*ny,dtype='int')]/16/10;
	azim0=data[np.linspace(1,nz*nx*ny*2-1,nz*nx*ny,dtype='int')]/10;

# 	dip0=dip0.reshape([nz,nx,ny],order='F')
# 	azim0=azim0.reshape([nz,nx,ny],order='F')
# 	plt.figure()
# 	plt.imshow(time0[:,:,0]);
# 	plt.jet()
# 	plt.colorbar();
# 	plt.title(name+'-'+str(time0.min())+'-'+str(time0.max()))
# 	plt.show()
	shot=np.array([0,0,shot[2]],dtype='float32');
	t,d,a=fmm.eikonal(vel,xyz=shot,ax=[ox,dx,nx],ay=[oy,dy,ny],az=[oz,dz,nz],order=1,angle=True);
	time=t.reshape(nx,ny,nz,order='F'); #[x,y,z]
	time=np.swapaxes(time,0,2); #[z,x,y]
# 	dip=d.reshape(nx,ny,nz,order='F');
# 	dip=np.swapaxes(dip,0,2); #[z,x,y]
# 	print("Angle difference: ",np.linalg.norm(dip0[:,:,0]-dip[:,0,:]))
	print("Time difference: ",np.linalg.norm(time0[:,:,0]-time[:,0,:]))
	fd=open(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm/layer.P.%s.time.buf'%name,'wb')
	np.float32(time).flatten(order='F').tofile(fd)
	np.float32(time).flatten(order='F').tofile(fd)
	
	
	fd=open(os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/timefmm/layer.P.%s.angle.buf'%name,'wb')
	data=np.zeros(nx*ny*nz*2,dtype=np.uint16)
	
	data[np.linspace(0,nz*nx*ny*2-2,nz*nx*ny,dtype='int')]=dip0*10*16*0;
	data[np.linspace(1,nz*nx*ny*2-1,nz*nx*ny,dtype='int')]=azim0*10*0;
	np.uint16(data).flatten(order='F').tofile(fd)
	
	
# 
# 
# 
# ## reference time from NLL
# file_time=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.AK_RC01_--.time.buf'
# # file_time=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.AK_CAPN_--.time.buf'
# 
# fd = open(file_time,'rb')
# time0=np.fromfile(fd, dtype = np.float32).reshape([nz,nx,ny*2],order='F')
# time0=time0[:,:,0].reshape([nz,nx,ny],order='F')
# 
# file_time=os.getenv('HOME')+'/DATALIB/softs/NonLinLoc/nlloc_sample_test/time/layer.P.AK_RC01_--.angle.buf'
# fd = open(file_time,'rb')
# data = np.fromfile(fd, dtype = np.uint16);#.reshape([nz,nx,ny*4],order='F')    #[zxy]
# dip0=data[np.linspace(0,nz*nx*ny*2-2,nz*nx*ny,dtype='int')]/16/10;
# azim0=data[np.linspace(1,nz*nx*ny*2-1,nz*nx*ny,dtype='int')]/10;
# 
# dip0=dip0.reshape([nz,nx,ny],order='F')
# azim0=azim0.reshape([nz,nx,ny],order='F')
# 

# 
# 
# ## pyekfmm
# vel=np.swapaxes(vel,0,2).reshape([nz*nx*ny,1],order='F'); #transpose to [xyz]
# shot=np.array([0,0,0],dtype='float32');
# t,d,a=fmm.eikonal(vel,xyz=shot,ax=[ox,dx,nx],ay=[oy,dy,ny],az=[oz,dz,nz],order=1,angle=True);
# time=t.reshape(nx,ny,nz,order='F'); #[x,y,z]
# time=np.swapaxes(time,0,2); #[z,x,y]
# 
# dip=d.reshape(nx,ny,nz,order='F');
# dip=np.swapaxes(dip,0,2); #[z,x,y]
# 
# azim=a.reshape(nx,ny,nz,order='F');
# azim=np.swapaxes(dip,0,2); #[z,x,y]
# 
# print("Angle difference: ",np.linalg.norm(dip0[:,:,0]-dip[:,0,:]))
# print("Time difference: ",np.linalg.norm(time0[:,:,0]-time[:,0,:]))
# 
# plt.figure();
# plt.subplot(1,3,1)
# plt.imshow(dip0[:,:,0],aspect='auto');
# plt.jet();
# plt.title('NLL')
# 
# 
# plt.subplot(1,3,2)
# plt.imshow(dip[:,0,:],aspect='auto');plt.gca().set_yticks([]);
# plt.jet();
# plt.title('pyekfmm')
# 
# plt.subplot(1,3,3)
# plt.imshow(dip0[:,:,0]-dip[:,0,:],clim=[0,50],aspect='auto');plt.gca().set_yticks([]);
# plt.jet();
# plt.title('NLL-pyekfmm')
# plt.colorbar(orientation='vertical',shrink=0.6);
# plt.show();
# 
# 
# 
# plt.figure();
# plt.subplot(1,3,1)
# plt.imshow(time0[:,:,0],aspect='auto');
# plt.jet();
# plt.title('NLL')
# plt.show()
# 
# plt.subplot(1,3,2)
# plt.imshow(time[:,0,:],aspect='auto');plt.gca().set_yticks([]);
# plt.jet();
# plt.title('pyekfmm')
# 
# plt.subplot(1,3,3)
# plt.imshow(time0[:,:,0]-time[:,0,:],clim=[0,200],aspect='auto');plt.gca().set_yticks([]);
# plt.jet();
# plt.title('NLL-pyekfmm')
# plt.colorbar(orientation='vertical',shrink=0.6);
# plt.show();












