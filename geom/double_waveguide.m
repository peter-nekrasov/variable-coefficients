function [coefs,dinds] = double_waveguide(X,Y,L,amp,len,wid,wdist,dsig,eps)

vpre = zeros(size(X));

ilow = (X >= -len/2) & (X <= len/2) & (Y >= -(wdist+wid)/2) & (Y <= -(wdist-wid)/2);
ihigh = (X >= -len/2) & (X <= len/2) & (Y <= (wdist+wid)/2) & (Y >= (wdist-wid)/2);

vpre(ilow) = 1;
vpre(ihigh) = 1;

N = size(X,1);
[src,targ,ind,sz,num] = get_fft_grid(N,L,0);

vpre_pad = [vpre zeros(N,N-1); zeros(N-1,2*N-1)];

kern = exp(-(targ(1,:).^2 + targ(2,:).^2)/(2*dsig.^2));
kern_aug = reshape(kern,size(vpre_pad));

ind = find(vecnorm(targ) < 1e-12);
[ii, jj] = ind2sub(size(kern_aug), ind);
kern_aug = circshift(kern_aug,-[ii-1,jj-1]);

kern_hat = fft2(kern_aug);
vpre_hat = fft2(vpre_pad);

V = ifft2(kern_hat.*vpre_hat);
V = V(1:N,1:N);

V = amp*V / max(V(:));

dinds = find(abs(V(:)) > eps );
coefs = V(dinds);

end