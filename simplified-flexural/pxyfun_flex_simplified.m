function [Kpxy,nbr] = pxyfun_flex_simplified(x,slf,nbr,l,ctr,quads,zk,V)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling helm_dirichlet_kernel

npxy = 100+round(zk*1.5*max(l)*2);
ths  = (0:(npxy-1))*(2*pi)/npxy;
rpxy = 1.5;
proxy = rpxy*[cos(ths);sin(ths)];

pxy = bsxfun(@plus,ctr,diag(l)*proxy);

rsrc = x(:,slf);
rtar = pxy;

Kpxy_helm = flexgreen(zk,rsrc,rtar);
Kpxy_helm = Kpxy_helm(:,:,1);

Kpxy_helm_1 = bsxfun(@times,Kpxy_helm,quads(slf).');
Kpxy_helm_trans_1 = bsxfun(@times,zk^2*V(slf).'/size(proxy,2),Kpxy_helm);

Kpxy = [Kpxy_helm_1;Kpxy_helm_trans_1];
ctruse = ctr(:);
dxy = abs(x(1:2,nbr)-ctruse)./l;
nbr = nbr(vecnorm(dxy) < 1.5);

end
