
clc;close all;clear;
addpath(genpath('/Users/chenyk/Downloads/toolbox_fast_marching/toolbox_fast_marching'));

fid=fopen('time_vgrad.bin','r');


d=fread(fid,[501,501],'float')';
t=d; 

figure;imagesc(t);colormap(jet);hold on;plot(1,1,'rp','markersize',20,'markerfacecolor','r');

z=[1:50:501];
end_points=[z;500*ones(size(z))];
paths = compute_geodesic(t,end_points);
for ii=1:length(paths)
plot(paths{ii}(2,:),paths{ii}(1,:),'g--','linewidth',2);
end




[tx,tz]=gradient(t);

path = stream3c([],[],[],tx,tz,[],end_points(1,1),end_points(2,1),[],0.1,10000)';







