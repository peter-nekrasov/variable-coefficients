function [val,grad,hess,third] = flexgreen(zk,src,targ)

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

r2 = (xt-xs).^2 + (yt-ys).^2;

if nargout == 1
val = hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;

val(r2 < 1e-8) = 1/(2*zk^2)*hkdiffgreen(zk,[0;0],[0;1e-8]);
elseif nargout == 2
[val,grad] = hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;
grad = 1/(2*zk^2)*grad;

val(r2 < 1e-8) = 1/(2*zk^2)*hkdiffgreen(zk,[0;0],[0;1e-8]);
elseif nargout == 3
[val,grad,hess] = hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;
grad = 1/(2*zk^2)*grad;
hess = 1/(2*zk^2)*hess;

val(r2 < 1e-8) = 1/(2*zk^2)*hkdiffgreen(zk,[0;0],[0;1e-8]);
elseif nargout == 4
[val,grad,hess,third] = hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;
grad = 1/(2*zk^2)*grad;
hess = 1/(2*zk^2)*hess;
third = 1/(2*zk^2)*third;

val(r2 < 1e-8) = 1/(2*zk^2)*hkdiffgreen(zk,[0;0],[0;1e-8]);
end

end