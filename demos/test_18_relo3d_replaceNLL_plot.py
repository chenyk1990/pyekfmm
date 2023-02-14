import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

file='/Users/chenyk/softs/GrowClust3D.jl/examples/data/out/out.ekfmmgrid3D.cat'
fid=open(file)
lines=fid.readlines()
print(len(lines))
#After Growclust
lons1=[ii.split()[8] for ii in lines]
lats1=[ii.split()[7] for ii in lines]
deps1=[ii.split()[9] for ii in lines]
lons1=np.array(lons1,dtype='float')
lats1=np.array(lats1,dtype='float')
deps1=np.array(deps1,dtype='float')
print('lons1',lons1.min(),lons1.max(),lons1.std(),lons1.mean())
print('lats1',lats1.min(),lats1.max(),lats1.std(),lats1.mean())
print('deps1',deps1.min(),deps1.max(),deps1.std(),deps1.mean())

#Before Growclust
lons2=[ii.split()[-2] for ii in lines]
lats2=[ii.split()[-3] for ii in lines]
deps2=[ii.split()[-1] for ii in lines]
lons2=np.array(lons2,dtype='float')
lats2=np.array(lats2,dtype='float')
deps2=np.array(deps2,dtype='float')
print('lons2',lons2.min(),lons2.max(),lons2.std(),lons2.mean())
print('lats2',lats2.min(),lats2.max(),lats2.std(),lats2.mean())
print('deps2',deps2.min(),deps2.max(),deps2.std(),deps2.mean())

## Using EKFMM
# lons1 -119.7194 -119.66167 0.004137240838628007 -119.69024502475246
# lats1 39.6405 39.6791 0.005028246631742764 39.66602459158416
# lons2 -119.7194 -119.66233 0.004326586284892841 -119.69024492574258
# lats2 39.6405 39.6805 0.005239158994711416 39.666024603960395
## Using NLL
# lons1 -119.7194 -119.66285 0.004118244857035998 -119.69024487623763
# lats1 39.6405 39.6791 0.005053245843940473 39.66602464108911
# lons2 -119.7194 -119.66233 0.004326586284892841 -119.69024492574258
# lats2 39.6405 39.6805 0.005239158994711416 39.666024603960395

## crosssection
# A
lon1=-119.698
lat1=39.658
# A'
lon2=-119.684
lat2=39.67

## plot
plt.figure(figsize=(16, 8))

plt.subplot(1,3,1);
plt.plot(lons2,lats2,'.',color='k');
plt.xlim([-119.700,-119.683]);
plt.ylim([39.676,39.653]);
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.locator_params(axis='x', nbins=6)
plt.title('Before relocation')
plt.xlabel('Lon (deg)');plt.ylabel('Lat (deg)');
plt.gca().invert_yaxis()
plt.plot([lon1,lon2],[lat1,lat2],'-',color='r');
plt.text(lon1-0.0005,lat1,'A',fontsize=14,color='r')
plt.text(lon2+0.0002,lat2,'A\' ',fontsize=14,color='r')
plt.text(-119.703,39.676, "a)", fontsize=24, color='k')

plt.subplot(1,3,2);
plt.plot(lons1,lats1,'.',color='k');plt.gca().invert_yaxis()
plt.xlim([-119.700,-119.683]);
plt.ylim([39.676,39.653]);
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.locator_params(axis='x', nbins=6)
plt.title('After relocation')
plt.xlabel('Lon (deg)');#plt.ylabel('Lat (deg)');
plt.gca().invert_yaxis()
plt.plot([lon1,lon2],[lat1,lat2],'-',color='r');
plt.text(lon1-0.0005,lat1,'A',fontsize=14,color='r')
plt.text(lon2+0.0002,lat2,'A\' ',fontsize=14,color='r')
plt.text(-119.703,39.676, "b)", fontsize=24, color='k')

# plt.show()

def point2line(x0,y0,x1,y1,x2,y2):
	k=(y2-y1)/(x2-x1);
	a=k;
	b=-1;
	c=y1-k*x1;
	
	x=(b*(b*x0-a*y0)-a*c)/(a*a+b*b);
	y=(a*(-b*x0+a*y0)-b*c)/(a*a+b*b);
	
	dist=np.sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y))
	return dist,(x,y)
	
##
#y=k(x-lon1)+lat1; 
#kx - y + lat1-k*lon1 = 0 -> ax+by+c=0;
#x=(b*(b*x0-a*y0)-a*c)/(a*a+b*b);
#y=(a*(-b*x0+a*y0)-b*c)/(a*a+b*b);
#
#https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#:~:text=In%20Euclidean%20geometry%2C%20the%20'distance,nearest%20point%20on%20the%20line.
#

##
deg2km=6371*2*np.pi/360.0
thr=0.1
xs1=[]
zs1=[]
for ii in range(len(lons1)):
	dist1=np.sqrt((lons1[ii]-lon1)*(lons1[ii]-lon1)+(lats1[ii]-lat1)*(lats1[ii]-lat1))
	dist2,tmp=point2line(lons1[ii],lats1[ii],lon1,lat1,lon2,lat2)
	dist2=dist2*deg2km;
	dist1=dist1*deg2km;
	if dist2<=thr:
		xs1.append(dist1)
		zs1.append(deps1[ii])
xs1=np.array(xs1)
zs1=np.array(zs1)

xs2=[]
zs2=[]
for ii in range(len(lons2)):
	dist1=np.sqrt((lons2[ii]-lon1)*(lons2[ii]-lon1)+(lats2[ii]-lat1)*(lats2[ii]-lat1))
	dist2,tmp=point2line(lons2[ii],lats2[ii],lon1,lat1,lon2,lat2)
	dist2=dist2*deg2km;
	dist1=dist1*deg2km;
	if dist2<=thr:
		xs2.append(dist1)
		zs2.append(deps2[ii])
xs2=np.array(xs2)
zs2=np.array(zs2)


plt.subplot(2,3,3);
plt.plot(xs2,zs2,'.',color='k');

plt.gca().set_xlim(xmin=0,xmax=1.8);
plt.gca().set_ylim(ymin=7.75,ymax=9.25);
#plt.xlabel('Distance (km)');
plt.ylabel('Depth (km)');
plt.title('Before relocation')
plt.gca().invert_yaxis();
plt.gca().text(-0.2,7.7, "c)", fontsize=24, color='k')

plt.subplot(2,3,6);
plt.plot(xs1,zs1,'.',color='k');
plt.gca().set_xlim(xmin=0,xmax=1.8);
plt.gca().set_ylim(ymin=7.75,ymax=9.25);
plt.xlabel('Distance (km)');plt.ylabel('Depth (km)');
plt.title('After relocation')
plt.gca().invert_yaxis();
plt.gca().text(-0.2,7.7, "d)", fontsize=24, color='k')

plt.savefig('exe_relo_ekfmm_new.png',format='png',dpi=300)
plt.savefig('exe_relo_ekfmm_new.pdf',format='pdf',dpi=300)

plt.show()
