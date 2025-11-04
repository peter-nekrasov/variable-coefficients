function val = flexgreen(zk,src,targ)

val = hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;

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

val(r2 < 1e-8) = 1/(2*zk^2)*hkdiffgreen(zk,[0;0],[0;1e-8]);

end