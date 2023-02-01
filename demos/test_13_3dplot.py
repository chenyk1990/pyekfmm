
## This DEMO is a 3D example [x,y,z] with constant velocity and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

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
    zdir='z', offset=Z.max(), **kw
)
_ = ax.contourf(
    X[0, :, :], data[0, :, :], Z[0, :, :],
    zdir='y', offset=0, **kw
)
C = ax.contourf(
    data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
    zdir='x', offset=X.max(), **kw
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
    xlabel='X [km]',
    ylabel='Y [km]',
    zlabel='Z [km]',
#     zticks=[0, -150, -300, -450],
)

# Set zoom and angle view
# ax.view_init(40, -30, 0)
# ax.set_box_aspect(None, zoom=0.9)

# Colorbar
fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Traveltime [s]')
plt.gca().scatter(0,0,0,marker='*r')

plt.savefig('test_13_plot3d.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)

# Show Figure
plt.show()

