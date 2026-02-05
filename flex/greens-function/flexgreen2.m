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

out = cat(3,val,hessxx,hessxy,hessyy,lap,gradlapx,gradlapy);

end