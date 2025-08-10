%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f = - k^2 V exp(i k x)
%
% Solved iteratively using FFT + GMRES 
%
%%%%%

L = 10; % length of grid
N1 = 300; % number of grid points

zk = 8;

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

[coefs,dinds] = bump(xxgrid,yygrid,-0.5,0.5,1,1e-8);
V = coefs(:,:,1);
coefs = coefs*zk^2;

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

% RHS (Incident field)
theta = -pi/3;
uinc = planewave(zk,[xxgrid(:) yygrid(:)].',theta);
rhs_vec = get_rhs(coefs,uinc,dinds);

figure(1); clf
tiledlayout(1,2);

nexttile
Vplot = zeros(length(xxgrid(:)),1);
Vplot(dinds) = V;
Vplot = reshape(Vplot,size(xxgrid));
pcolor(xxgrid,yygrid,Vplot); shading interp;
colorbar
title('V')
drawnow

nexttile
rhsplot = zeros(length(xxgrid(:)),1);
rhsplot(dinds) = rhs_vec;
rhsplot = reshape(rhsplot,size(xxgrid));
pcolor(xxgrid,yygrid,real(rhsplot)); shading interp;
colorbar
title('rhs')
drawnow

% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
[inds,corrs] = get_correct_helm(h,zk);
gfunc = @(s,t) helmgreen1(zk,s,t);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
% kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);

% Solve with GMRES
start = tic;
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,spmats,h,dinds,iinds,jinds,N2),rhs_vec,[],1e-10,200);
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('Time to solve: %5.2e s\n',t1)

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

% usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
usca = usca(:,:,1);

utot = usca + uinc;
utot = reshape(utot,size(xxgrid));

figure(2);
tiledlayout(1,3)

nexttile
pc = pcolor(xxgrid,yygrid,real(mu)); shading interp;
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
err = get_fin_diff_err(xxgrid,yygrid,utot,h,coefs,0.1,0.1,zk,dinds,'helm');

fprintf('Finite difference error: %.4e \n',err)

return