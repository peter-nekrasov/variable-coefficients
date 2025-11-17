function [Kpxy,nbr] = pxyfun_flex_simplified(x,slf,nbr,l,ctr,quads,zk,V)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling helm_dirichlet_kernel

npxy = 100+round(zk*1.5*max(l)*2);
ths  = (0:(npxy-1))*(2*pi)/npxy;
proxy = [1.5*cos(ths) 1.75*cos(ths) ; ...
    1.5*sin(ths) 1.75*sin(ths) ];

nx = [cos(ths) cos(ths)].'; 
ny = [sin(ths) sin(ths)].';

pxy = bsxfun(@plus,ctr,diag(l)*proxy);

rsrc = x(:,slf);
rtar = pxy;

[Kpxy_flex, Kpxy_flex_grad] = flexgreen(zk,rsrc,rtar);
Kpxy_flex_grad = nx.*Kpxy_flex_grad(:,:,1) + ny.*Kpxy_flex_grad(:,:,2);

Kpxy_flex_1 = bsxfun(@times,Kpxy_flex,quads(slf).');
Kpxy_flex_grad_1 = bsxfun(@times,Kpxy_flex_grad,quads(slf).');
Kpxy_helm_trans_1 = bsxfun(@times,-zk^4*V(slf).'/size(proxy,2),Kpxy_flex);
Kpxy_helm_grad_trans_1 = bsxfun(@times,-zk^4*V(slf).'/size(proxy,2),Kpxy_flex_grad);

Kpxy = [Kpxy_flex_1;Kpxy_flex_grad_1;Kpxy_helm_trans_1;Kpxy_helm_grad_trans_1];
ctruse = ctr(:);
dxy = abs(x(1:2,nbr)-ctruse)./l;
nbr = nbr(vecnorm(dxy) < 1.75);

end
