function [Kpxy,nbr] = pxyfun_fg(x,slf,nbr,l,ctr,quads,zk,rts,ejs,coefs)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling dirichlet_kernel

npxy = 100+round(zk*1.5*max(l)*2);
rmin = 1.5;
rmax = rmin+max(2*pi/zk/max(l),0.5*rmin);
proxy = get_proxy_ann(rmin,5,npxy,rmax);

pxy = bsxfun(@plus,ctr,diag(l)*proxy);

rsrc = x(:,slf);
rtar = pxy;

Kpxy_helm = fggreen(rsrc,rtar,rts,ejs);

Kpxy_1 = Kpxy_helm(:,:,1);
Kpxy_1 = bsxfun(@times,Kpxy_1,quads(slf).');

dfac = vecnorm(pxy - ctr(:)).^3*(1/rmin-1/rmax)/(max(l))/10;
Kpxy_trans_1 = bsxfun(@times, permute(coefs(slf,:,:),[2 1 3]),Kpxy_helm);
Kpxy_trans_1 = sum(Kpxy_trans_1,3);
Kpxy_trans_1 = bsxfun(@times, dfac.', Kpxy_trans_1);

Kpxy = [Kpxy_1;Kpxy_trans_1];
ctruse = ctr(:);
dxy = abs(x(1:2,nbr)-ctruse)./l;
nbr = nbr(vecnorm(dxy) < 1.5);

end
