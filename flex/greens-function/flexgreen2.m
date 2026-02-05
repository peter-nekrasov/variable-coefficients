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

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dy2 = dy.*dy;
r2 = dx2 + dy2;

val(r2<1e-16) = 0;
hessxx(r2<1e-16) = 0;
hessxy(r2<1e-16) = 0;
hessyy(r2<1e-16) = 0;
lap(r2<1e-16) = 0;
gradlapx(r2<1e-16) = 0;
gradlapy(r2<1e-16) = 0;

out = cat(3,val,hessxx,hessxy,hessyy,lap,gradlapx,gradlapy);

end