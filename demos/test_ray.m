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

