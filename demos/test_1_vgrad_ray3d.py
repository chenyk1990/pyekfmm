
## This DEMO is a 2D example [x,z] with velocity gradient and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np
import matplotlib.pyplot as plt

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');
t=fmm.eikonal(vel,xyz=np.array([0.5,0,0]),ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2);
time=t.reshape(101,101,101,order='F'); #[x,y,z]


## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',0.49965078, 0.0, 0.08905013, 0.007929926])


import matplotlib.pyplot as plt
import numpy as np

# Define dimensions
Nx, Ny, Nz = 101, 101, 101
X, Y, Z = np.meshgrid(np.arange(Nx)*0.01, np.arange(Ny)*0.01, np.arange(Nz)*0.01)

# Specify the 3D data
data=np.transpose(time,(1,0,2)); ## data requires [y,x,z] so tranpose the first and second axis

kw = {
    'vmin': data.min(),
    'vmax': data.max(),
    'levels': np.linspace(data.min(), data.max(), 10),
}

# Create a figure with 3D ax
fig = plt.figure(figsize=(8, 8))
plt.jet()
ax = fig.add_subplot(111, projection='3d')

# Plot contour surfaces
_ = ax.contourf(
    X[:, :, -1], Y[:, :, -1], data[:, :, -1],
    zdir='z', offset=Z.max(), alpha=0.5, **kw
)
_ = ax.contourf(
    X[0, :, :], data[0, :, :], Z[0, :, :],
    zdir='y', offset=0, alpha=0.5, **kw
)
C = ax.contourf(
    data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
    zdir='x', offset=X.max(), alpha=0.5, **kw
)
# --


# Set limits of the plot from coord limits
xmin, xmax = X.min(), X.max()
ymin, ymax = Y.min(), Y.max()
zmin, zmax = Z.min(), Z.max()
ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

# Plot edges
# edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
# ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
# ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
# ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

# Set labels and zticks
ax.set(
    xlabel='X (km)',
    ylabel='Y (km)',
    zlabel='Z (km)',
#     zticks=[0, -150, -300, -450],
)

# Set zoom and angle view
# ax.view_init(40, -30, 0)
# ax.set_box_aspect(None, zoom=0.9)

# Colorbar
fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, format= "%4.2f", label='Traveltime (s)')

plt.gca().scatter(0.5,0,0,s=500,marker='*',color='r')
plt.gca().set_xlim(0,1);
plt.gca().set_ylim(0,1);
plt.gca().set_zlim(0,1);

plt.savefig('test_1_vgrad_ray3d.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)

# Show Figure
plt.show()




# 
# 
# v1=1;
# v2=3;
# nz=501;
# nx=501;
# dx=0.01;
# dz=0.01;
# # vel=3.0*np.ones([501*501,1],dtype='float32'); #velocity axis must be x,y,z respectively
# v=np.linspace(v1,v2,nz);
# v=np.expand_dims(v,1);
# h=np.ones([1,nx])
# vel=np.multiply(v,h,dtype='float32'); #z,x
# # plt.imshow(vel);plt.jet();plt.show()
# 
# t=fmm.eikonal(vel.transpose().flatten(order='F'),xyz=np.array([0,0,0]),ax=[0,dx,nx],ay=[0,0.01,1],az=[0,dz,nz],order=2);
# time=t.reshape(nx,nz,order='F');#first axis (vertical) is x, second is z
# time=time.transpose(); #z,x
# 
# # tz=np.gradient(time,axis=1);
# # tx=np.gradient(time,axis=0);
# # # or
# tz,tx = np.gradient(time)
# 
# 
# 
# 
# plt.figure();
# plt.imshow(time,cmap=plt.cm.jet, interpolation='none', extent=[0,5,5,0]); #transpose so that first axis is z, second is x
# plt.plot(0,0.0,'*r',markersize=20);
# plt.xlabel('X (km)');plt.ylabel('Z (km)');
# plt.jet()
# plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
# 
# for ii in range(1,502,50):
# 	paths,nrays=fmm.stream2d(-tx,-tz, 501, ii, step=0.1, maxvert=10000)
# 	plt.plot(500*dx,(ii-1)*dz,'vb',markersize=15);
# 	## plot rays
# 	plt.plot((paths[0,:]-1)*dx,(paths[1,:]-1)*dz,'g--',markersize=20);
# 
# plt.savefig('test_1_vgrad_ray.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.savefig('test_1_vgrad_ray.pdf',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
# 
# plt.show()
# 
# 
# 
# 
# 
# ## Verify
# print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
# print(['Correct result:',4.4039335, 0.0, 0.7773146, 0.604218])
# 


