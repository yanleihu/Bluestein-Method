function g=Bluestein_dft(x,f1,f2,fs,mout)

[m,n]=size(x);
f11=f1+(mout*fs+f2-f1)/(2*mout);
f22=f2+(mout*fs+f2-f1)/(2*mout);
a = exp(j*2*pi*f11/fs);
w = exp(-j*2*pi*(f22-f11)/(mout*fs)); 
h=((-m+1):max(mout-1,m-1)).';
np=(0:(m-1))';
wp=w.^((h.^2)./2);
ap=a.^(-np).*wp(m+np);
xx=x.*ap(:,ones(1,n));
nf=2^nextpow2(m+mout-1);
ftt=fft(1./wp(1:(m+mout-1)),nf);
ft=fft(xx,nf);
ft=ft.*ftt(:,ones(1,n));
g=ifft(ft);
g=g(m:(m+mout-1),:).*wp(m:(m+mout-1),ones(1,n));

l=linspace(0,mout-1,mout);l=l./mout.*(f22-f11)+f11;
Mshift=-m/2; 
Mshift=repmat(exp(-1i.*2*pi.*l.*(Mshift+1/2)/fs),[n 1]);
g=g.'.*Mshift;

end
