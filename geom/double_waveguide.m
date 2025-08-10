function coefs = double_waveguide(X,Y,amp,len,wid,wdist,dsig)

vpre = zeros(size(X));

ilow = (X >= -len/2) & (X <= len/2) & (Y >= -(wdist+wid)/2) & (Y <= -(wdist-wid)/2);
ihigh = (X >= -len/2) & (X <= len/2) & (Y <= (wdist+wid)/2) & (Y >= (wdist-wid)/2);

vpre(ilow) = 1;
vpre(ihigh) = 1;

N = size(X,1);

vpre_pad = [vpre zeros(N); zeros(N,2*N)];

kern = exp(-(X.^2 + Y.^2)/(2*dsig.^2));
kern_pad = [kern zeros(N); zeros(N,2*N)];

[ii,jj] = find((abs(X) < 1e-12) & (abs(Y) < 1e-12));
kern_pad = circshift(kern_pad,-[ii,jj]);

kern_hat = fft2(kern_pad);
vpre_hat = fft2(vpre_pad);

V = ifft2(kern_hat.*vpre_hat);
V = V(1:N,1:N);

V = amp*V / max(V(:));

coefs = {V};

end