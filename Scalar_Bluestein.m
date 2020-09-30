%%%%%%Scalar diffraction computation method using Bluestein method
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

function [gout,pixelout] = Scalar_Bluestein(gin,mxin,myin,pixelin,d,xstart,xend,ystart,yend,mxout,myout)
%%%%%% Defination of input variable:
%%%%%% gin----Complex amplitude of the incident light beam
%%%%%% mxin,myin-----resolution of the input plane in transverse and longitudinal directions
%%%%%% pixelin----pixel size of the incident light beam
%%%%%% d----diffraction distance
%%%%%% xstart,xend----computation range in transverse direction
%%%%%% ystart,yend----computation range in longitudinal direction
%%%%%% mxout,myout----resolution of the output plane in transverse and longitudinal directions
%%%%%% Defination of output variable:
%%%%%% gout----Complex amplitude of the outgoing light beam
%%%%%% pixelout----pixel size of the outgoing light beam

global lamda k

L0=(myin-1)*pixelin;                                                        % dimension of the diffraction plane
x0=linspace(-L0/2,L0/2,mxin);                                               % transverse coordinates (x) in the diffraction plane
y0=linspace(-L0/2,L0/2,myin);                                               % longitudinal coordinates (y) in the diffraction plane
[x0,y0]=meshgrid(x0,y0);

pixelout=(xend-xstart)/(mxout-1);                                           % pixel size of the outgoing light beam
x1=linspace(xstart,xend,mxout);                                             % transverse coordinates (x) in the imaging plane
y1=linspace(ystart,yend,myout);                                             % longitudinal coordinates (y) in the imaging plane
[x1,y1]=meshgrid(x1,y1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          % calculating scalar diffraction below
F0=exp(j*k*d)/(j*lamda*d).*exp(j*k/2/d*(x1.^2+y1.^2));                      % assign exp(ikd)/(i¦Ëd)exp[ik(x2+y2) /2d]; Equation 4 in the paper
F=exp(j*k/2/d*(x0.^2+y0.^2));                                               % assign exp[ik (x02+y02) /2d]; Equation 5 in the paper
gout=gin.*F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       % using Bluestein method to calculate the complex amplitude of the outgoing light beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            % one-dimensional FFT in one direction
fs = lamda*d/pixelin;                                                       % dimension of the imaging plane
fy1=ystart+fs/2;
fy2=yend+fs/2;
gout = Bluestein_dft(gout,fy1,fy2,fs,myout);%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            % one-dimensional FFT in the other direction
fx1=xstart+fs/2;
fx2=xend+fs/2;
gout = Bluestein_dft(gout,fx1,fx2,fs,mxout);%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gout=F0.*gout;                                                              % obtain the complex amplitude of the outgoing light beam
end
