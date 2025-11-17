%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% nonlinear simplified flexural scattering problem:
%
%      \Delta^2 u - k^4 (1 + V(x) (a0 + a1 | u |^2 ) ) u = 0
%
%
% Solved iteratively using poor man's Newton
%
%%%%%

% incfield = 'pointsource';
incfield = 'planewave';

a0 = 1;
a1 = 0.5;
zk = 2;

step_size = 0.01;

L = 20; % length of grid
N1 = 250; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

[coefs,dinds] = bump(xxgrid,yygrid,1.5,1,1,1e-8);
V = coefs(:,:,1);

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

% RHS (Incident field)

if strcmpi(incfield,'planewave')

theta = -pi/4;
uinc = planewave(zk,[xxgrid(:) yygrid(:)].',theta);

elseif strcmpi(incfield,'pointsource')

uinc = flexgreen(zk,[-L/3;-L/3],[xxgrid(:) yygrid(:)].');

end

rhs_pde = zk.^4*coefs.*(a0 + a1*uinc(dinds,:,:).*conj(uinc(dinds,:,:))).*uinc(dinds,:,:);

figure(1); clf
tiledlayout(1,2);

nexttile
Vtot = zeros(length(xxgrid(:)),1);
Vtot(dinds) = V;
Vtot = reshape(Vtot,size(xxgrid));
pcolor(xxgrid,yygrid,Vtot); shading interp;
colorbar
title('V')
drawnow

nexttile
rhsplot = zeros(length(xxgrid(:)),1);
rhsplot(dinds) = rhs_pde;
rhsplot = reshape(rhsplot,size(xxgrid));
pcolor(xxgrid,yygrid,real(rhsplot)); shading interp;
colorbar
title('rhs')
drawnow

% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,0);
[inds,corrs] = get_correct_flex_simplified(h);
gfunc = @(s,t) flexgreen(zk,s,t);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
% idspmats = id_plus_corr_sum(coefs,spmats,dinds,h);
% kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);

mu_k = zeros(size(dinds));
utot_k = zeros(size(dinds));

as = 0:step_size:a1;

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

% rhs_pde = zk.^4*coefs.*(a0 + a1*uinc(dinds,:,:).*conj(uinc(dinds,:,:))).*uinc(dinds,:,:);

for ii = 1:numel(as)

disp(ii)

a1_k = as(ii)

ie_resid = 1;

while vecnorm(ie_resid) > 1e-8

rhs_vec = zk^4*V.*(a0 + a1_k*(utot_k).*conj(utot_k)).*uinc(dinds,:,:)  ;
rhsnorm = vecnorm(rhs_vec);

coefs = -zk.^4*V.*(a0+a1_k.*utot_k.*conj(utot_k));

% Solve with GMRES
start = tic;
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,idspmats,iinds,jinds,N2),rhs_vec,[],1e-10,200);
mu_k = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-12,5000,[],[],mu_k);
t1 = toc(start);

% usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
u_k = sol_eval_fft(mu_k,evalkerns,iinds,jinds,N1,N2);
u_k = u_k(dinds,:,1);
utot_k = u_k + uinc(dinds,:,:);

ie_resid = mu_k - zk^4*V.*(a0 + a1_k*(utot_k).*conj(utot_k)).*utot_k ;

end

end

%%

mu_tot = zeros(size(xxgrid));
mu_tot(dinds) = mu_k;

usca = sol_eval_fft(mu_k,evalkerns,iinds,jinds,N1,N2);

utot = usca + uinc;
utot = reshape(utot,size(xxgrid));

figure(2);
tiledlayout(1,3)

nexttile
pc = pcolor(xxgrid,yygrid,real(mu_tot)); shading interp;
title('Re(\mu)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(utot)); shading interp;
title('Re(u)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot)); shading interp;
title('|u|')
colorbar
       
% Calculate error with finite difference

pdecoef = -zk^4*Vtot.*(a0 + a1*utot.*conj(utot));
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,pdecoef,0.1,0.1,zk,dinds,'flex-nonlinear');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)