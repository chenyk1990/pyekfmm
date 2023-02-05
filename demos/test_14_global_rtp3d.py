## This DEMO is a 3D spherical example [r,t,p] with constant velocity and with one shot
#  in a global scale 
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt

##Event lat/lon/dep
evlat=-25+90
evlon=0+180
evdep=0

## size of rtp axes (number of samples)
numr=1
numt=100*2
nump=200*2

## beg/end in r
begr=6400
endr=6400

begt=0
endt=180

## beg/end in p
begp=0
endp=360

## sampling in r
dr=(endr-begr)/numr
dt=(endt-begt)/numt
dp=(endp-begp)/nump


vel=6.0*np.ones([numr*numt*nump,1],dtype='float32');
t=fmm.eikonal_rtp(vel,rtp=np.array([6400-evdep,evlat,evlon]),ar=[begr,dr,numr],at=[begt,dt,numt],ap=[begp,dp,nump],order=1);#spherical needs order=1
time=t.reshape(numr,numt,nump,order='F'); #[r,t,p]

# import matplotlib.pyplot as plt
# # plt.figure(figsize=(12, 7))
# plt.imshow(time[numr-1,:,:],extent=[begp-180,endp-180,endt-90,begt-90],aspect='auto')#,clim=(0, 2000.0)); #rtp->tp
# plt.gca().invert_yaxis()
# plt.xlabel('Lon (deg)'); 
# plt.ylabel('Lat (deg)');
# plt.jet()
# plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
# plt.plot(evlon-180,evlat-90,'*',color='r',markersize=15)
# # plt.gca().invert_yaxis()
# plt.savefig('test_14_global.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.savefig('test_14_global.pdf',format='pdf',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.show()
# 
# ## Verify
# print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
# print(['Correct result:',2140.937, 0.016649984, 75.647644, 5722.5664])
# #print(['Correct result:',285.1792, 0.016649984, 75.70713, 5731.57]) #If using order=2
# 

fig = plt.figure(figsize=(8, 8))
plt.jet()
ax = fig.add_subplot(projection='3d')

# Make data
u = np.linspace(0, 2 * np.pi, nump)
v = np.linspace(0, np.pi, numt)

x = 10 * np.outer(np.cos(u), np.sin(v))
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

from matplotlib import cm
my_col = cm.jet(np.squeeze(time.transpose())/np.amax(time))

# Plot the surface
cp=ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors = my_col, alpha=0.2)
plt.axis('off')

ax.plot(0, 0, 0,'*',color='r',markersize=15)

# Set an equal aspect ratio
# ax.set_aspect('equal')
cb = plt.colorbar(cp, cax = fig.add_axes([0.37,0.1,0.3,0.02]), format= "%4.0f", orientation='horizontal',label='Traveltime (s)',values=np.linspace(0,np.amax(time),200))
plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# X = np.arange(-5, 5, 0.25)
# Y = np.arange(-5, 5, 0.25)
# X, Y = np.meshgrid(X, Y)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)
# my_col = cm.jet(np.random.rand(Z.shape[0],Z.shape[1]))
# 
# surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors = my_col,
#         linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)
