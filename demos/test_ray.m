clc;close all;clear;

fid=fopen('time.bin','r');


d=fread(fid,[501,501],'float')';
figure;imagesc(d);colormap(jet);hold on;plot(251,1,'rp','markersize',20,'markerfacecolor','r');
t=d;

s=[1,251];

r=[375,375];

plot(r(2),r(1),'bv','markersize',20,'markerfacecolor','b');

dz=1;dx=1;
rr=r; %current
pp=r; %previous
dray=sqrt(dx^2+dz^2)*2;
while norm(s-rr)>=2
% for ii=1:10
    ddz=(t(rr(1)+1,rr(2))-t(rr(1),rr(2)))/dz;
    ddx=(t(rr(1),rr(2)+1)-t(rr(1),rr(2)))/dx;
%     tant=-ddz/ddx;
    theta=atan(-ddz/ddx);
%     theta/pi*180
%     round(dray*cos(theta)*sign(theta))
%     round(dray*sin(theta))
    pp=rr;
%     rr=[rr(1)+round(dx*tant/dz),rr(2)+sign(tant)];
    rr=[rr(1)+round(dray*sin(theta)),rr(2)+round(dray*cos(theta)*sign(theta))];
    plot([pp(2) rr(2)],[pp(1) rr(1)],'g--','linewidth',2);
end
% end



%% 
clc;close all;clear;
fid=fopen('/Users/chenyk/Downloads/tomo_codes/vsmooth5063.bin','r');
vel=fread(fid,[50,63],'float');
figure;imagesc(vel);colormap(jet);colorbar;hold on;

fid=fopen('timev2d.bin','r');
d=fread(fid,[63,50],'float')';
% figure;imagesc(d);colormap(jet);hold on;
plot(11,11,'rp','markersize',20,'markerfacecolor','r');
t=d;

s=[11,11];
dz=1;dx=1;
dray=sqrt(dx^2+dz^2)*2;
for rz=15:4:40
r=[rz,40];
rr=r; %current
while norm(s-rr)>=1.5
% for ii=1:200
    ddz=0.5*((t(rr(1)+1,rr(2))-t(rr(1),rr(2)))/dz + (t(rr(1),rr(2))-t(rr(1)-1,rr(2)))/dz);
    ddx=0.5*((t(rr(1),rr(2)+1)-t(rr(1),rr(2)))/dx + (t(rr(1),rr(2))-t(rr(1),rr(2)-1))/dx);
%     tant=-ddz/ddx;
    theta=abs(atan(ddz/ddx));
%     theta/pi*180
%     round(dray*cos(theta)*sign(theta))
%     round(dray*sin(theta))
    pp=rr;
%     rr=[rr(1)+round(dx*tant/dz),rr(2)+sign(tant)];
rr1=rr(1)-round(dray*sin(theta)*sign(ddz));
rr2=rr(2)-round(dray*cos(theta)*sign(ddx));
    
%     if rr1<=1
%         rr1=1;
%     end
%     if rr1>500;
%         rr1=500;
%     end
%     if rr2<=1
%         rr2=1;
%     end
%     if rr2>500;
%         rr2=500;
%     end

    rr=[rr1,rr2];
    plot([pp(2) rr(2)],[pp(1) rr(1)],'g--','linewidth',2);
end

end


figure;imagesc(vel);colormap(jet);hold on;plot(s(2),s(1),'rp','markersize',20,'markerfacecolor','r');

z=[15:4:40];
end_points=[z;40*ones(size(z))];
paths = compute_geodesic(t,end_points);
for ii=1:length(paths)
plot(paths{ii}(2,:),paths{ii}(1,:),'g--','linewidth',2);
end



figure;imagesc(vel);colormap(jet);hold on;plot(s(2),s(1),'rp','markersize',20,'markerfacecolor','r');
[tx,tz]=gradient(t);
for ii=[15:4:40];
path2 = yc_stream3c(-tx,-tz,40.5,ii,0.1,10000);

path2 = yc_trimrays(path2,[11,11],0.2);

plot(path2(1,:),path2(2,:),'g--','linewidth',2);
end














