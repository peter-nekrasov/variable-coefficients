function [src,targ,ind,sz,num] = get_fft_grid(N,L,itwon)

if nargin < 3
    itwon = 0;
end

if itwon 

M = 2^(nextpow2(N));
xl = L*(-M:M-1)/(N-1);

else

% xl = L*(-(N-1):(N-1))/(N-1);
xl = linspace(-L,L,2*N-1);

end

[XL,YL] = meshgrid(xl);

src = [0;0];
targ = [XL(:).'; YL(:).'];

ind = find((XL == 0) & (YL ==0));
sz = size(XL);

num = length(xl);

end