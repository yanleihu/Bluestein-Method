%%%%%%Scalar calculation of the converging spherical wave
%%%%%%The work is published by Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%Please cite properly if you use the codes for any use
%%%%%%unit: um

clear all;
clc;
tic;

global lamda k
lamda=800e-3;
k=2*pi/lamda;

my0=1081;
mx0=1081;
pixel0=8;
L0=(mx0-1)*pixel0;
[xx,yy]=meshgrid(-(my0-1)/2:(my0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
Aperture=sign(1-sign(xx.^2+yy.^2-((my0-1)./2).^2));
A0=Aperture;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=600e3;%
g=A0.*exp(-1i.*k./2./f.*(x0.^2+y0.^2));
% figure
% imshow(angle(g),[]);
L=0.2e3;
ly=201;
zstart=571e3;
zfinish=629e3;
num=round((zfinish-zstart)/lamda+1);
num=round(num/ly);

zstart=600e3-(ly-1).*num.*lamda./2;
zfinish=600e3+(ly-1).*num.*lamda./2;

x1start=0;
x1end=0;
y1start=-L./2;
y1end=L./2; 
mx1=1;
my1=1081;
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);                   
% h=waitbar(0,'Crosection Caculating ...');
index=1;
for d=linspace(zstart,zfinish,ly)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I=abs(g1).^2;
    P=angle(g1);
    I_z(:,index)=I(:,round((mx1+1)/2));
    P_z(:,index)=P(:,round((mx1+1)/2));
    index=index+1;
%     waitbar((index-1)./ly);
end

% close(h);
I_z=I_z./max(I_z(:));
[z1,y1]=meshgrid(zstart:((zfinish-zstart)./(ly-1)):zfinish,linspace(y1start,y1end,my1));
toc
figure
surf(z1./1e6,y1./1e3,I_z), colormap hot,axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-0.1:0.1:0.1],'FontSize',24)  

figure
surf(z1./1e6,y1./1e3,P_z), colormap hot,axis tight, view([0, -270]), colorbar('Ticks',[-3,0,3]), shading interp;
set(gca,'YTick',[-0.1:0.1:0.1],'FontSize',24)   
