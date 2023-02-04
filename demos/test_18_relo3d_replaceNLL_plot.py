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
lons1=np.array(lons1,dtype='float')
lats1=np.array(lats1,dtype='float')
print('lons1',lons1.min(),lons1.max(),lons1.std(),lons1.mean())
print('lats1',lats1.min(),lats1.max(),lats1.std(),lats1.mean())

#Before Growclust
lons2=[ii.split()[-2] for ii in lines]
lats2=[ii.split()[-3] for ii in lines]
lons2=np.array(lons2,dtype='float')
lats2=np.array(lats2,dtype='float')
print('lons2',lons2.min(),lons2.max(),lons2.std(),lons2.mean())
print('lats2',lats2.min(),lats2.max(),lats2.std(),lats2.mean())

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

## plot
plt.figure(figsize=(16, 8))

plt.subplot(1,2,1);
plt.plot(lons2,lats2,'.',color='k');
plt.xlim([-119.700,-119.683]);
plt.ylim([39.676,39.653]);
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.title('Before relocation')
plt.xlabel('Lon (deg)');plt.ylabel('Lat (deg)');
plt.gca().invert_yaxis()

plt.subplot(1,2,2);
plt.plot(lons1,lats1,'.',color='k');plt.gca().invert_yaxis()
plt.xlim([-119.700,-119.683]);
plt.ylim([39.676,39.653]);
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.title('After relocation')
plt.xlabel('Lon (deg)');plt.ylabel('Lat (deg)');
plt.gca().invert_yaxis()
plt.savefig('exe_relo_ekfmm.png',format='png',dpi=300)
plt.savefig('exe_relo_ekfmm.pdf',format='pdf',dpi=300)
plt.show()


