import numpy as np
import matplotlib.pyplot as plt

fid=open('/Users/chenyk/softs/NonLinLoc/nlloc_sample_test/loc-old/alaska.sum.grid0.loc.hyp','r')
lines=fid.readlines()

ne=0;
for ii in range(len(lines)):
	line=lines[ii].split()
	if len(line) != 0 and line[0]=='STAT_GEOG' and float(line[6])>0:
		print(lines[ii])
		ne=ne+1;
	
lons=np.zeros(ne);
lats=np.zeros(ne);
deps=np.zeros(ne);

ne=0;
for ii in range(len(lines)):
	line=lines[ii].split()
	if len(line) != 0 and line[0]=='STAT_GEOG' and float(line[6])>0:
		lons[ne]=float(line[4])
		lats[ne]=float(line[2])
		deps[ne]=float(line[6])
		ne=ne+1;
		




fid=open('/Users/chenyk/softs/NonLinLoc/nlloc_sample_test/locfmm/alaska.sum.grid0.loc.hyp','r')
lines=fid.readlines()

ne=0;
for ii in range(len(lines)):
	line=lines[ii].split()
	if len(line) != 0 and line[0]=='STAT_GEOG' and float(line[6])>0:
		print(lines[ii])
		ne=ne+1;
	
lons2=np.zeros(ne);
lats2=np.zeros(ne);
deps2=np.zeros(ne);

ne=0;
for ii in range(len(lines)):
	line=lines[ii].split()
	if len(line) != 0 and line[0]=='STAT_GEOG' and float(line[6])>0:
		lons2[ne]=float(line[4])
		lats2[ne]=float(line[2])
		deps2[ne]=float(line[6])
		ne=ne+1;
		


fid=open('/Users/chenyk/softs/NonLinLoc/nlloc_sample_test/obs/station_coordinates.txt');
lines=fid.readlines();
lines=lines[1:]
stnames=[ii.split()[1] for ii in lines];

stlons=[float(ii.split()[4]) for ii in lines];
stlats=[float(ii.split()[3]) for ii in lines];
stdeps=[float(ii.split()[5]) for ii in lines];

lonmin=-156;lonmax=-143;
latmin=59;latmax=64;
lonmin=-152;lonmax=-148;
latmin=61;latmax=62;

stlons2=stlons
stlats2=stlats
stdeps2=stdeps
stlons=[]
stlats=[]
stdeps=[]
for ii in range(len(stlons2)):
	if stlons2[ii]>lonmin and stlons2[ii]<lonmax and stlats2[ii]>latmin and stlats2[ii]<latmax:
		stlons.append(stlons2[ii])
		stlats.append(stlats2[ii])
		stdeps.append(stdeps2[ii])

# eids2=[ii for ii in eids if ii not in eids0]#there are 703 events in eids,

fig=plt.figure(figsize=(8,8))
cax = fig.add_axes([0.1,0.3,0.6,0.6])
plt.plot(lons,lats,'r*',alpha=0.7);
plt.plot(lons2,lats2,'g*',alpha=0.7);	
plt.plot(stlons,stlats,'kv');	
plt.gca().set_xlim(xmin=lonmin,xmax=lonmax);
plt.gca().set_ylim(ymin=latmin,ymax=latmax);	
plt.setp(plt.gca().get_xticklabels(), visible=False)
plt.ylabel('Latitude (deg)')
plt.plot(-149.8997,61.2176,'*y',label='Manual',markersize=14);
plt.text(-149.8997,61.2176,'Anchorage, Alaska',color='y',fontsize=14)


#below is to add frameboxes to highlight the comparison
import matplotlib.patches as patches
for ii in range(len(lons)):
	rect = patches.Rectangle((lons[ii]-0.05, lats[ii]-0.0125), 0.1, 0.025, linewidth=1, edgecolor='b', facecolor='none')
	plt.gca().add_patch(rect)

cax = fig.add_axes([0.1,0.1,0.6,0.15])
plt.plot(lons,deps,'r*',alpha=0.7);	
plt.plot(lons2,deps2,'g*',alpha=0.7);	
plt.plot(stlons,stdeps,'kv');	
plt.gca().invert_yaxis()
plt.gca().set_xlim(xmin=lonmin,xmax=lonmax);
plt.xlabel('Longitude (deg)')
plt.ylabel('Depth (km)')
#below is to add frameboxes to highlight the comparison
for ii in range(len(lons)):
	rect = patches.Rectangle((lons[ii]-0.05, deps[ii]-3), 0.1, 6, linewidth=1, edgecolor='b', facecolor='none')
	plt.gca().add_patch(rect)
	
	
cax = fig.add_axes([0.75,0.3,0.15,0.6])
plt.plot(deps,lats,'r*',alpha=0.7);	
plt.plot(deps2,lats2,'g*',alpha=0.7);	
plt.plot(stdeps,stlats,'kv');
plt.setp(plt.gca().get_yticklabels(), visible=False)
plt.gca().set_ylim(ymin=latmin,ymax=latmax);	
plt.xlabel('Depth (km)')
#below is to add frameboxes to highlight the comparison
for ii in range(len(lons)):
	rect = patches.Rectangle((deps[ii]-3, lats[ii]-0.0125), 6, 0.025, linewidth=1, edgecolor='b', facecolor='none')
	plt.gca().add_patch(rect)


plt.legend(['NLL', 'FMM','Station'],loc='best', bbox_to_anchor=(0.8, -0.35, 0.2, 0.2))
plt.savefig('test_11_comp2NLL_2d_angle_location_plot.png',format='png',dpi=300)
plt.savefig('test_11_comp2NLL_2d_angle_location_plot.pdf',format='pdf',dpi=300)
plt.show()