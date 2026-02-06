% solve the forward problem given G and corrections

%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural waves
%
%%%%%

L = 1000; % length of grid
N1 = 300; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

[~,~,pcoefs] = bump4(xxgrid,yygrid,0.5,50,4,1e-12);

a0 = pcoefs{1}; 
b0 = pcoefs{3}; 
zk = (b0/a0)^0.25;
clear pcoefs

% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
[inds,corrs] = get_correct_flex(h);
gfunc = @(s,t) flexgreen2(zk,s,t);
%kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);
kerns = kerns / a0;

amp = 0.4;
[coefs0,dinds0,pcoefs,alpha0] = bump4(xxgrid,yygrid,amp,50,4,1e-12);
[iinds,jinds] = ind2sub(size(xxgrid),dinds0);



fprintf('Number of points: %d \n',size(dinds0,1))

% RHS (Incident field)

incfield = 'pointsource';

uincs = flexgreen2(zk,[-L/2;-L/2],[xxgrid(:) yygrid(:)].');
uinc = uincs(:,:,1);


rhs_vec0 = get_rhs(coefs0,uincs(dinds0,:,:));

start = tic;
[usca0, sol0] = solve_flex(rhs_vec0, coefs0, kerns, iinds, jinds, N1, N2);
mu = zeros(size(xxgrid));
mu(dinds0) = sol0;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

coefspad = zeros(numel(xxgrid),1,size(coefs0,3));
coefspad(dinds0,:,:) = coefs0;


deltas = 10.^(-3:-1:-4);
dus = zeros(2,length(deltas));
for i = 1:length(deltas)
delta = deltas(i);
[coefs,dinds,~,alpha] = bump4(xxgrid,yygrid,amp+delta,50,4,1e-12);
rhs_vec = get_rhs(coefs,uincs(dinds,:,:));
[iinds,jinds] = ind2sub(size(xxgrid),dinds);
[usca, sol] = solve_flex(rhs_vec, coefs, kerns, iinds, jinds, N1, N2);

alphascale = alpha0(dinds)./alpha(dinds);

% dcoefs = coefs - reshape(coefspad(dinds,:),[],1,size(coefs,3)).*alphascale;


[coefsp,dindsp,~,alphap] = bump4(xxgrid,yygrid,amp+delta,50,4,1e-12);
[coefsm,dindsm,~,alpham] = bump4(xxgrid,yygrid,amp-delta,50,4,1e-12);
% dcoefs = (coefsp.*alphap(dinds) - coefsm.*alpham(dinds))./alpha(dinds)/2/delta;
coefsmpad = zeros(numel(alpham(:)),1,size(coefs,3));
coefsmpad(dindsm,:,:) = coefsm.*alpham(dindsm);
dcoefs = (coefsp.*alphap(dindsp) - coefsmpad(dindsp,:,:))/2/delta;
dalpha = (alphap-alpham)/2/delta;


dcoefs = (coefsp - coefsm)/2/delta;

% coefsmpad(dindsm,:,:) = coefsm;
% dcoefs = (coefsp - coefsmpad(dindsp,:,:))/2/delta;
rhs2 = fast_apply_fft(mu(dinds),kerns,dcoefs,iinds,jinds,N2);
rhs2 = rhs2 - mu(dinds);

rhs2 = usca(dinds)./alpha(dinds);

% rhs2 = rhs2 + (dalpha(dindsp)).*mu(dinds);
% rhs2 = rhs2./alpha(dinds);
[freshu, ~] = solve_flex(-rhs2, coefs, kerns, iinds, jinds, N1, N2);
% [freshu2, ~] = solve_flex(-(dalpha(dindsp)).*mu(dinds)./alpha(dinds), coefs, kerns, iinds, jinds, N1, N2);
dus(1,i) = norm(usca - usca0) / norm(usca);
dus(2,i) = norm(usca - freshu*delta- usca0) / norm(usca);

% [norm(coefsp.*alphap(dindsp)  - coefspad(dindsp,:,:).*alpha0(dindsp),'fro')/norm(coefsp.*alphap(dindsp),'fro')...
% ,norm(coefsp.*alphap(dindsp) - delta*dcoefs  - coefspad(dindsp,:,:).*alpha0(dindsp),'fro')/norm(coefsp.*alphap(dindsp),'fro'),...
% norm(alphap - alpha0), norm(alphap - delta*dalpha - alpha0)]

% norm(dcoefs,'fro')
end

[deltas; dus]

utot = uinc + usca0;

utot = reshape(utot,size(xxgrid));

figure(2);
pc = pcolor(xxgrid,yygrid,real(mu)); shading interp;
title('Re(\mu)')
colorbar

figure(3); clf
tiledlayout(1,2)

nexttile
pc = pcolor(xxgrid,yygrid,real(utot)); shading interp;
title('real(u)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot)); shading interp;
title('|u|')
colorbar
       
% Calculate error with finite difference
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,pcoefs,10,10,zk,dinds0,'flex');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)

return



function [usca, sol] = solve_flex(rhs_vec, coefs, kerns, iinds, jinds, N1, N2)
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,200);
evalkerns = kerns(:,:,1);
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
end