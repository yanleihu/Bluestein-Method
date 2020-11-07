%Calculation of the optical fields using Debye theory.
%unit: um

clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 201;
global lamda k n1 NA fo
n1 = 1.514;
z = 0e0;
Mo = 100;
F = 180e3;
N=(M-1)/2; 
lamda = 800e-3;
NA = 1.4; 
fo=F./Mo;
k = (2.*pi)./lamda;
polar='rc';
m=linspace(-M/2,M/2,M);
n=linspace(-M/2,M/2,M);
[m n] = meshgrid(m,n); 
th=asin(NA.*sqrt(m.^2 + n.^2)./(N.*n1));
thh=th;
th(thh>asin(NA./n1))=0;
phi = atan2 (n,m);
phi(phi<0) = phi(phi<0)+2.*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fi=mod(phi.*1,2*pi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0 = ones(M);
E = E0.*exp(i.*Fi);
E(thh>asin(NA./n1))=0;
figure
imshow(angle(E),[]); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my=501;
mx=501;
startx = -0.8e0; 
endx = 0.8e0;  
starty = -0.8e0; 
endy = 0.8e0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ex Ey Ez] = Vector_Bluestein(E,M,polar,startx,endx,starty,endy,z,mx,my); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ix = abs(Ex).^2;
Iy = abs(Ey).^2;
Iz = abs(Ez).^2;
Phasex=angle(Ex);
Phasey=angle(Ey);
Phasez=angle(Ez);
I = Ix+Iy+Iz;
Imax=max(max(I));
I=I./Imax;
Ix=Ix./Imax;
Iy=Iy./Imax;
Iz=Iz./Imax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x = linspace(startx,endx,mx);
y = linspace(starty,endy,my);
figure
surfc(x./lamda,y./lamda,I), colormap hot, axis equal, view([0, -270]),colorbar,shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24) 
% saveas(gcf,'cxyIntensity.png');
figure
surfc(x./lamda,y./lamda,Ix), colormap hot, axis equal, view([0, -270]),colorbar,shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24)   
% saveas(gcf,'cxyIntensityx.png');
figure
surfc(x./lamda,y./lamda,Iy), colormap hot, axis equal, view([0, -270]),colorbar,shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24)   
% saveas(gcf,'cxyIntensityy.png');
figure
surfc(x./lamda,y./lamda,Iz), colormap hot, axis equal, view([0, -270]),colorbar,shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24) 
% saveas(gcf,'cxyIntensitz.png');
figure
surfc(x./lamda,y./lamda,Phasex), colormap hot, axis equal, view([0, -270]),colorbar('Ticks',[-3,0,3]),shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24)   
% saveas(gcf,'cxyphasex.png');
figure
surfc(x./lamda,y./lamda,Phasey), colormap hot, axis equal, view([0, -270]),colorbar('Ticks',[-3,0,3]),shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24)  
% saveas(gcf,'cxyphasey.png');
figure
surfc(x./lamda,y./lamda,Phasez), colormap hot, axis equal, view([0, -270]),colorbar('Ticks',[-3,0,3]),shading interp;
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24)   
% saveas(gcf,'cxyphasez.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startz=-0.8e0;
endz=0.8e0;
mz=100;
my=100;
mx=1;
II=zeros(my,mz);
h=waitbar(0,'Caculating XZ plane...');
index=0;
loop=mz;
for z=linspace(startz,endz,mz)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [Ex Ey Ez] = Vector_Bluestein(E,M,polar,0,0,starty,endy,z,mx,my); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ix = abs(Ex).^2;
Iy = abs(Ey).^2;
Iz = abs(Ez).^2;
I = Ix+Iy+Iz;
index=index+1;
waitbar(index./loop)
II(:,index)=I(:,round((mx+1)/2));
end
close(h);
IImax=max(max(II));

z=linspace(startz,endz,mz);
y=linspace(starty,endy,my);
figure
surfc(z./lamda,y./lamda,II./IImax), colormap hot, axis equal, view([0, -270]),colorbar,
shading interp
set(gca,'YTick',[-1:1:1],'XTick',[-1:1:1],'FontSize',24)
% saveas(gcf,'cxzIntensity.png');
