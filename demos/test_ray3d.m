clc;close all;clear;

fid=fopen('time3d.bin','r');


d=fread(fid,[101,101*101],'float')';
d=reshape(d,101,101,101);%original: x,y,z
d=yc_transp(d,23);
% d=yc_transp(d,13);

% figure;imagesc(d);colormap(jet);hold on;

frames=[1,101,1];
yc_mada3d(d,frames,1:101,1:101,1:101);colormap(jet);caxis([0,0.5]);alpha(0.5);

hold on;
plot3(51,1,1,'rp','markersize',20,'markerfacecolor','r');

s=[1,51,1];%z,x,y
r=[75,75,75];
plot3(r(1),r(2),r(3),'bv','markersize',20,'markerfacecolor','b');


t=d;
dz=1;dx=1;dy=1;
rr=r; %current,z,x,y
pp=r; %previous,z,x,y
dray=sqrt(dx^2+dz^2+dy^2)*1;
while norm(s-rr)>=1
% for ii=1:1
    ddz=(t(rr(1)+1,rr(2),rr(3))-t(rr(1),rr(2),rr(3)))/dz;
    ddx=(t(rr(1),rr(2)+1,rr(3))-t(rr(1),rr(2),rr(3)))/dx;
    ddy=(t(rr(1),rr(2),rr(3)+1)-t(rr(1),rr(2),rr(3)))/dy;
%     tant=-ddz/ddx;

    ddr=sqrt(ddz^2+ddx^2+ddy^2);
    ddr2=sqrt(ddx^2+ddy^2);
    theta=acos(ddr2/ddr);
    phi=atan(-ddy/ddx);
    theta/pi*180
    phi/pi*180
%     theta=atan(-ddz/ddx);
%     theta/pi*180
%     round(dray*cos(theta)*sign(theta))
%     round(dray*sin(theta))
    pp=rr;
%   rr=[rr(1)+round(dray*sin(theta)),rr(2)+round(dray*cos(theta)*sign(theta))];

    rr1=rr(1)-round(dray*sin(theta));
    rr2=rr(2)+round(dray*cos(theta)*cos(phi)*sign(phi));
    rr3=rr(3)+round(dray*cos(theta)*sin(phi));

    if rr1<=1
        rr1=1;
    end
    if rr1>101;
        rr1=101;
    end
    if rr2<=1
        rr2=1;
    end
    if rr2>101;
        rr2=101;
    end
    if rr3<=1
        rr3=1;
    end
    if rr3>101;
        rr3=101;
    end
    rr=[rr1,rr2,rr3];
    
    plot3([pp(2) rr(2)],[pp(3) rr(3)],[pp(1) rr(1)],'g--','linewidth',2);
end


%%
addpath(genpath('/Users/chenyk/Downloads/toolbox_fast_marching/toolbox_fast_marching'));

figure;
frames=[1,101,1];
yc_mada3d(d,frames,1:101,1:101,1:101);colormap(jet);caxis([0,0.5]);alpha(0.5);

hold on;
plot3(51,1,1,'rp','markersize',20,'markerfacecolor','r');

s=[1,51,1];%z,x,y
r=[75,10.2,75];
plot3(r(1),r(2),r(3),'bv','markersize',20,'markerfacecolor','b');


z=[1:50:501];
end_points=[r(3);r(1);r(2)]; %z,x,y in endpoints
paths = compute_geodesic(t,end_points);%z,x,y in paths | z,x,y in t
plot3(paths(2,:),paths(3,:),paths(1,:),'g--','linewidth',2);


% for ii=1:length(paths)
% plot(paths{ii}(3,:),paths{ii}(2,:),paths{ii}(1,:),'g--','linewidth',2);
% end




