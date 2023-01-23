
## This DEMO is a 3D example [x,y,z] with constant velocity and with one shot
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin

import pyekfmm as fmm
import numpy as np

vel=3.0*np.ones([101*101*101,1],dtype='float32');
t=fmm.eikonal(vel,xyz=np.array([0.5,0,0]),ax=[0,0.01,101],ay=[0,0.01,101],az=[0,0.01,101],order=2);
time=t.reshape(101,101,101,order='F'); #[x,y,z]


import matplotlib.pyplot as plt
plt.imshow(time[:,50,:].transpose()); #xyz->xz, transpose so that first axis is z, second is x
plt.xlabel('X');plt.ylabel('Z');
plt.jet()
plt.colorbar(orientation='horizontal',shrink=0.6);
plt.show()


## Verify
print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
print(['Correct result:',0.49965078, 0.0, 0.08905013, 0.007929926])



import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def plot_3D_array_slices(array):
    min_val = array.min()
    max_val = array.max()
    n_x, n_y, n_z = array.shape
    colormap = plt.cm.YlOrRd

    x_cut = array[n_x//2,:,:]
    Y, Z = np.mgrid[0:n_y, 0:n_z]
    X = n_x//2 * np.ones((n_y, n_z))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2,facecolors=colormap((x_cut-min_val)/(max_val-min_val)), shade=False)
    ax.set_title("x slice")

    y_cut = array[:,n_y//2,:]
    X, Z = np.mgrid[0:n_x, 0:n_z]
    Y = n_y//2 * np.ones((n_x, n_z))
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2, facecolors=colormap((y_cut-min_val)/(max_val-min_val)), shade=False)
    ax.set_title("y slice")

    z_cut = array[:,:,n_z//2]
    X, Y = np.mgrid[0:n_x, 0:n_y]
    Z = n_z//2 * np.ones((n_x, n_y))
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2, facecolors=colormap((z_cut-min_val)/(max_val-min_val)), shade=False)
    ax.set_title("z slice")

    plt.show()


n_pts = 100
r_square = (np.mgrid[-1:1:1j*n_pts, -1:1:1j*n_pts, -1:1:1j*n_pts]**2).sum(0)
plot_3D_array_slices(r_square)