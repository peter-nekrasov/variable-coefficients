function out = flexgreen2(zk,src,targ)

[val,~,hess,third] = hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;
% grad = 1/(2*zk^2)*grad;
hess = 1/(2*zk^2)*hess;
third = 1/(2*zk^2)*third;

hessxx = hess(:,:,1);
hessxy = hess(:,:,2);
hessyy = hess(:,:,3);

lap = hessxx+hessyy;
gradlapx = third(:,:,1) + third(:,:,3);
gradlapy = third(:,:,2) + third(:,:,4);

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dy2 = dy.*dy;
r2 = dx2 + dy2;

eulgam = 0.5772156649015328606065120900824;

val(r2<1e-16) = 1/(2*zk^2)/(2*pi)*(log(2/zk)-log(2/(1i*zk)));
hessxx(r2<1e-16) = 1/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
hessxy(r2<1e-16) = 0;
hessyy(r2<1e-16) = 1/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
lap(r2<1e-16) = 2/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
gradlapx(r2<1e-16) = 0;
gradlapy(r2<1e-16) = 0;

out = cat(3,val,hessxx,hessxy,hessyy,lap,gradlapx,gradlapy);

end