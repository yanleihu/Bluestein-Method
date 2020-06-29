%%%Full path calculation of light propagation
%%%
%%%%unit: um

clear all;
clc;
tic;

global lamda k n1 NA fo
lamda=800e-3;                                                               % wavelength
k=2*pi/lamda;                                                               % propagation (wave) constant free space
n1 = 1.514;                                                                 % refractive index of oil-immersion objective lens
NA = 1.4;                                                                   % numerical aperture of the objective
Mo = 100;                                                                   % magnification of the lens
F = 180e3;                                                                  % tube length
fo=F./Mo;                                                                   % objective focal length (F=tube Length, M=magnification)
R=fo.*NA./n1;                                                               % objective back aperture (radius)

my0=255;                                                                   
mx0=255;                                                                   % resolution of the input plane
pixel0=8;                                                                   % pixel size
L0=(mx0-1)*pixel0;                                                          % dimension of input plane
[xx,yy]=meshgrid(-(my0-1)/2:(my0-1)/2,-(my0-1)/2:(my0-1)/2);
Aperture=sign(1-sign(xx.^2+yy.^2-((my0-1)./2).^2));                         % circular aperture
A0=Aperture;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                             % here define the input phase profile, e.g. a hologram or a continuous phase 
% g=double(imread('CGHEE255.bmp'));g=g./255.*2.*pi;                        % hologram 
g=3.*atan2(yy,xx);
g=A0.*exp(1i.*g);                                                           % 3-fold helical vortex
figure
imshow(angle(g),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate the light after propagation d
d=600e3;                                                                    % propagation distance
L=4*L0;                                                                     % desired dimension of imaging plane
x1start=-L./2;                                                              % start positon x
x1end=L./2;                                                                 % end position x
y1start=-L./2;                                                              % start positon y
y1end=L./2;                                                                 % end position y
mx1=255;                                                                   % desired resolution of imaging plane
my1=255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=linspace(x1start,x1end,mx1); 
y1=linspace(y1start,y1end,my1);
[x1,y1]=meshgrid(x1,y1); 
% figure
% surf(x1,y1,abs(g1).^2), axis equal, axis tight, view([0, 270]), colorbar, shading interp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% propagation through a lens
f=600e3;                                                                    % lens focal lenght 600mm
g1=g1.*exp(-1i.*k./2./f.*(x1.^2+y1.^2));                                    % transmittance function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate the light after propagation d
d=800e3;
L=2.*L0;
x2start=-L./2;
x2end=L./2;
y2start=-L./2;
y2end=L./2;
mx2=255;
my2=255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel2]=Scalar_Bluestein(g1,mx1,my1,pixel1,d,x2start,x2end,y2start,y2end,mx2,my2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=linspace(x2start,x2end,mx2); 
y2=linspace(y2start,y2end,my2); 
[x2,y2]=meshgrid(x2,y2); 
% figure
% surf(x2,y2,abs(g1).^2), axis equal, axis tight, view([0, 270]), colorbar, shading interp;
% figure
% imshow(angle(g1),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%propagation through a lens
f=200e3;                                                                    % lens focal length
g1=g1.*exp(-1i.*k./2./f.*(x2.^2+y2.^2));                                    % transmittance function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%alculate the light after propagation d
d=200e3+1.8e3;
L=2*R;
x3start=-L./2;
x3end=L./2;
y3start=-L./2;
y3end=L./2;
mx3=255;
my3=255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel3]=Scalar_Bluestein(g1,mx2,my2,pixel2,d,x3start,x3end,y3start,y3end,mx3,my3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3=linspace(x3start,x3end,mx3); 
y3=linspace(y3start,y3end,my3);
[x3,y3]=meshgrid(x3,y3);
figure
surf(x3,y3,abs(g1).^2), axis equal, axis tight, view([0, 270]), colorbar, shading interp;
figure
imshow(angle(g1),[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%vector diffraction xy plane
M = 255;                                                                   % resolution of input plane
my=255;                                                                     % desired resolution in the imaging plane
mx=255;
startx = -5*lamda;                                                          % start position x
endx = 5*lamda;                                                             % end position x
starty = -5*lamda;                                                          % start position y
endy = 5*lamda;                                                             % end position y
z = 0e0;                                                                    % position shift along the optical axis
%%%%%%%%%%%%%%                                                              % polarization; corresponds to apodization functions, here Sine condition is used;
% polar='x';
% polar='y';
polar='rc';
% polar='lc';
% polar='ra';
% polar='az';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           % input light amplitude and phase
E=g1;
% E(thh>asin(NA./n1))=0;                                                    % remove parts outside numerical aperture
% figure
% surf(m, n,angle(E)), title 'on the obj. pup.',axis equal, axis tight, view([0, 270]), colorbar, shading interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        % vector diffraction using Bluestein method
[Ex Ey Ez] = Vector_Bluestein(E,M,polar,startx,endx,starty,endy,z,mx,my);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%                                                            % calculation of intensity for each component
Ix = abs(Ex).^2;
Iy = abs(Ey).^2;
Iz = abs(Ez).^2;
I = Ix+Iy+Iz;                                                               % total intensity
Imax=max(max(I));
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                % calibrate field dimensions
x = linspace(startx,endx,mx);
y = linspace(starty,endy,my);
figure
surfc(x./lamda,y./lamda,I./Imax), colormap hot, axis equal, view([180, -90]),shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   % vector diffraction yz plane
startz=-40*lamda;                                                            % start position z
endz=40*lamda;                                                               % end position z
startx =0;
endx = 0;
starty = -20*lamda;                                                          % start position y
endy = 20*lamda;                                                             % end position y
mz=50;                                                                      % desired resolution z
my=100;                                                                     % desired resolution y
mx=1;

II=zeros(my,mz);
h=waitbar(0,'Caculating XZ plane...');
index=0;
loop=mz;
for z=linspace(startz,endz,mz)
    [Ex Ey Ez] = Vector_Bluestein(E,M,polar,startx,endx,starty,endy,z,mx,my);
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
surfc(z./lamda,y./lamda,II./IImax), colormap hot, axis tight, view([180, -90]),
shading interp

toc
