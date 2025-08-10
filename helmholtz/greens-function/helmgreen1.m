function kerns = helmgreen1(zk,src,targ)
%CHNK.HELM2D.GREEN evaluate the helmholtz green's function
% for the given sources and targets

[~,ns] = size(src);
[~,nt] = size(targ);

eulgam = 0.5772156649015328606065120900824;

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

rx = xt-xs;
ry = yt-ys;

rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);

kerns = cell(1);

if nargout > 0
    h0 = besselh(0,1,zk*r);
    val = 0.25*1i*h0;

    val(r < 1e-12) = 1/(2*pi)*(1i*pi/2  - eulgam + log(2/zk));
end

kerns{1} = -val;
