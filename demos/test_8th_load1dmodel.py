## This DEMO is a 3D spherical example [r,t,p] with input 1D velocity and with one shot
#  in a regional scale (e.g., Texas)
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

f=open("data/ex_mymodel.nd");#example velocity model, compatible to pytaup 1D model file
lines=f.readlines()[0:25]; #cut before mantle, Moho depth is 43 km.
zs=[float(line.strip().split(' ')[0]) for line in lines]
vps=[float(line.strip().split(' ')[1]) for line in lines]

# y=np.interp(np.linspace(0,zmax,numr), zs, vps)

##Event lat/lon/dep
evla=54.344124
evlo=-117.248193
evdp=3.270

##Target station lat/lon/dep
stla=54.3256
stlo=-117.201103
stel=0

## size of rtp axes (number of samples)
numr=201
numt=201
nump=201

## beg/end in r (radius/depth)
zmax=3.270
endr=6371
begr=endr-zmax

## beg/end in t (latitude)
begt=54.3256
endt=54.344124   #ev

## beg/end in p (longitude)
begp=-117.248193 #ev
endp=-117.201103

## sampling in r
dr=(endr-begr)/numr
dt=(endt-begt)/numt
dp=(endp-begp)/nump


## station index [numr,0,nump]
irst = numr-1
itst = 0
ipst = nump-1

## construct 3D velocity model
vp1d=np.flipud(np.interp(np.linspace(0,zmax,numr), zs, vps)); #from zmax to zmin/from rmin to rmax
vel=6*np.ones([numr,numt,nump],dtype='float32');
for ir in range(numr):
	for it in range(numt):
		for ip in range(nump):
			vel[ir,it,ip]=vp1d[ir];
vel=vel.reshape([numr*numt*nump,1],order='F');
# vel=3*np.ones([numr*numt*nump,1],dtype='float32');
## construct 3D velocity model

t=fmm.eikonal_rtp(vel,rtp=np.array([endr-evdp,evla,evlo]),ar=[begr,dr,numr],at=[begt,dt,numt],ap=[begp,dp,nump],order=2);#spherical needs order=1
time=t.reshape(numr,numt,nump,order='F'); #[r,t,p]
print('P arrival time:',time[irst,itst,ipst])
print('Correct arrival time:',1.5274389)


import matplotlib.pyplot as plt
plt.imshow(time[numr-1,:,:],extent=[round(begp*100)/100,round(endp*100)/100,endt,begt],aspect='auto')#,clim=(0, 2000.0)); #rtp->tp
plt.gca().invert_yaxis()
plt.xlabel('Lon'); 
plt.ylabel('Lat');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.plot(evlo,evla,'*',color='r',markersize=15)
plt.plot(stlo,stla,'v',color='b',markersize=15)
plt.show()


##Cartesian is correct
vp1d=np.interp(np.linspace(0,zmax,numr), zs, vps); #from zmax to zmin/from rmin to rmax
vel=6*np.ones([numt,nump,numr],dtype='float32');
for ir in range(numr):
	for it in range(numt):
		for ip in range(nump):
			vel[it,ip,ir]=vp1d[ir];
vel=vel.reshape([numr*numt*nump,1],order='F');
# vel=3*np.ones([numr*numt*nump,1],dtype='float32');
tc=fmm.eikonal(vel,xyz=np.array([evla*111.1949,evlo*111.1949,evdp]),ax=[begt*111.1949,dt*111.1949,numt],ay=[begp*111.1949,dp*111.1949,nump],az=[0,dr,numr],order=2);#spherical needs order=1
timec=tc.reshape(numt,nump,numr,order='F'); #[x,y,z]

print('Cartesian P arrival time:',timec[0,nump-1,0])
print('Correct P arrival time:',1.5274389)



# """
# Example of takeoff angle calculation using obspy taup module
# """
# from obspy.taup import TauPyModel
# from obspy.taup import taup_create
# from obspy.geodetics.base import locations2degrees
# 
# 
# taup_create.build_taup_model('./data/ex_mymodel.nd',output_folder='./data/')
# model = TauPyModel(model="./data/ex_mymodel.npz")
# # calculate travel time
# dist=locations2degrees(evla,evlo,stla,stlo)
# arrivals=model.get_travel_times(source_depth_in_km=evdp+stel,distance_in_degree=dist,
#                                 phase_list=("p"),receiver_depth_in_km=0)
# for arr in arrivals:
#     print(arr.ray_param, arr.time, arr.takeoff_angle)
#     
# # calculate ray path
# arrivals=model.get_ray_paths(source_depth_in_km=evdp+stel,distance_in_degree=dist,
#                                 phase_list=("p"),receiver_depth_in_km=0)
# arrivals.plot_rays(plot_type="cartesian")
# 
# 
# ## Verify
# print('P arrival time:',time[irst,itst,ipst])
# # print(['P arrival time from taup:',time[irst,itst,ipst])
# 



