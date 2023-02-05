## This DEMO is a 3D spherical example [r,t,p] with constant velocity and with one shot
#  in a global scale 
# 
#  COPYRIGHT: Yangkang Chen, 2022, The University of Texas at Austin


import matplotlib.pyplot as plt
import numpy as np

#...some data processing
used_theta=np.linspace(0,np.pi,100)
used_rad=np.linspace(6200,6400,100)
data2D=np.random.randn(100,100)
theta,rad = np.meshgrid(used_theta, used_rad) #rectangular plot of polar data
X = theta
Y = rad
# 
# fig = plt.figure()
# ax = fig.add_subplot(111,projection='polar')
# ax.pcolormesh(X, Y, data2D) #X,Y & data2D must all be same dimensions
# 
# ax.set_thetamin(45)
# ax.set_thetamax(135)
# # ax.set_radmin(6200)
# 
# ax.set_rorigin(6000)
# # ax.set_theta_zero_location('W', offset=10)
# 
# plt.show()
# 


import pyekfmm as fmm


##Event lat/lon/dep
evlat=-25+90
evlon=0+180
evdep=400

## size of rtp axes (number of samples)
numr=2000
numt=1
nump=2000

## beg/end in r
begr=3200
endr=6400

begt=evlat
endt=evlat

## beg/end in p
begp=90
endp=270

## sampling in r
dr=(endr-begr)/numr
dt=(endt-begt)/numt
dp=(endp-begp)/nump


vel=6.0*np.ones([numr*numt*nump,1],dtype='float32');
t=fmm.eikonal_rtp(vel,rtp=np.array([6400-evdep,evlat,evlon]),ar=[begr,dr,numr],at=[begt,dt,numt],ap=[begp,dp,nump],order=1);#spherical needs order=1
time=t.reshape(numr,nump,order='F'); #[r,t,p]


## plot on polar coordinates
used_rad=np.linspace(begr,endr,numr)
used_theta=np.linspace((begp-180)/180*np.pi,(endp-180)/180*np.pi,nump)
theta,rad = np.meshgrid(used_theta, used_rad) #rectangular plot of polar data
X = theta
Y = rad

fig = plt.figure()
plt.jet()
ax = fig.add_subplot(111,projection='polar')
cm=ax.pcolormesh(X, Y, time) #X,Y & data2D must all be same dimensions

ax.plot(0, 6400-evdep,'*',color='r',markersize=15)

ax.set_thetamin(begp)
ax.set_thetamax(endp)
ax.set_thetamin(-90)
ax.set_thetamax(+90)
ax.set_rorigin(5800-5700)
ax.set_theta_zero_location('W', offset=90)
ax.set_theta_zero_location('N', offset=0)
# ax.set_rmax(2)
ax.set_rticks([3200])  # less radial ticks
ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
ax.grid(True)
# ax.set_xlabel('Longitude (deg)')
ax.set_ylabel('Longitude (deg)')
# ---- mod here ---- #
# ax.set_theta_zero_location("N")  # theta=0 at the top
ax.set_theta_direction(-1)  # theta increasing clockwise

cb = plt.colorbar(cm, cax = fig.add_axes([0.37,0.1,0.3,0.02]), format= "%4.0f", orientation='horizontal',label='Traveltime (s)')


plt.savefig('test_14_global_rp.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.savefig('test_14_global_rp.pdf',format='pdf',dpi=300,bbox_inches='tight', pad_inches=0)

plt.show()

# fig = plt.figure()
# plt.imshow(time)
# plt.jet()
# plt.show()





# import matplotlib.pyplot as plt
# plt.imshow(time[numr-1,:,:],extent=[begp-180,endp-180,endt-90,begt-90],aspect='auto')#,clim=(0, 2000.0)); #rtp->tp
# plt.gca().invert_yaxis()
# plt.xlabel('Lon (deg)'); 
# plt.ylabel('Lat (deg)');
# plt.jet()
# plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
# plt.plot(evlon-180,evlat-90,'*',color='r',markersize=15)
# # plt.gca().invert_yaxis()
# plt.savefig('test_14_global_rp.png',format='png',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.savefig('test_14_global_rp.pdf',format='pdf',dpi=300,bbox_inches='tight', pad_inches=0)
# plt.show()

## Verify
# print(['Testing result:',time.max(),time.min(),time.std(),time.var()])
# print(['Correct result:',2140.937, 0.016649984, 75.647644, 5722.5664])
#print(['Correct result:',285.1792, 0.016649984, 75.70713, 5731.57]) #If using order=2



