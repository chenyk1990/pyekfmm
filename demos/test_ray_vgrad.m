clc;close all;clear;
addpath(genpath('/Users/chenyk/Downloads/toolbox_fast_marching/toolbox_fast_marching'));

fid=fopen('time_vgrad.bin','r');


d=fread(fid,[501,501],'float')';
t=d; 
figure;imagesc(t);colormap(jet);hold on;plot(1,1,'rp','markersize',20,'markerfacecolor','r');

s=[1,1];

r=[500,500];%

plot(r(2),r(1),'bv','markersize',10,'markerfacecolor','b');

dz=1;dx=1;
rr=r; %current
pp=r; %previous
dray=sqrt(dx^2+dz^2)*2;
while norm(s-rr)>=2
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
% end

%% ray 2
% r=[450,500];
% rr=r; %current
% while norm(s-rr)>=2
% % for ii=1:10
%     ddz=(t(rr(1)+1,rr(2))-t(rr(1),rr(2)))/dz;
%     ddx=(t(rr(1),rr(2)+1)-t(rr(1),rr(2)))/dx;
% %     tant=-ddz/ddx;
%     theta=atan(-ddz/ddx);
% %     theta/pi*180
% %     round(dray*cos(theta)*sign(theta))
% %     round(dray*sin(theta))
%     pp=rr;
% %     rr=[rr(1)+round(dx*tant/dz),rr(2)+sign(tant)];
%     rr=[rr(1)+round(dray*sin(theta)),rr(2)+round(dray*cos(theta)*sign(theta))];
%     plot([pp(2) rr(2)],[pp(1) rr(1)],'g--','linewidth',2);
% end

% 
% r=[300,500];
% rr=r; %current
% while norm(s-rr)>=2
% % for ii=1:10
%     ddz=(t(rr(1)+1,rr(2))-t(rr(1),rr(2)))/dz;
%     ddx=(t(rr(1),rr(2)+1)-t(rr(1),rr(2)))/dx;
% %     tant=-ddz/ddx;
%     theta=atan(-ddz/ddx);
% %     theta/pi*180
% %     round(dray*cos(theta)*sign(theta))
% %     round(dray*sin(theta))
%     pp=rr;
% %     rr=[rr(1)+round(dx*tant/dz),rr(2)+sign(tant)];
%     rr=[rr(1)+round(dray*sin(theta)),rr(2)+round(dray*cos(theta)*sign(theta))];
%     plot([pp(2) rr(2)],[pp(1) rr(1)],'g--','linewidth',2);
% end
% 
%% ray N
for rz=2:50:450
r=[rz,500];
rr=r; %current
while norm(s-rr)>=2
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

% print(gcf,'-dpng','-r300','test1.png');

figure;imagesc(t);colormap(jet);hold on;plot(1,1,'rp','markersize',20,'markerfacecolor','r');

z=[1:50:501];
end_points=[z;500*ones(size(z))];
paths = compute_geodesic(t,end_points);
for ii=1:length(paths)
plot(paths{ii}(2,:),paths{ii}(1,:),'g--','linewidth',2);
end


%% my version
[tx,tz]=gradient(t);

figure;imagesc(t);colormap(jet);hold on;plot(1,1,'rp','markersize',20,'markerfacecolor','r');
for ii=1:50:501
path2 = yc_stream3c(-tx,-tz,500,ii,10,1000)';

plot(path2(:,1),path2(:,2),'g--','linewidth',2);
end







