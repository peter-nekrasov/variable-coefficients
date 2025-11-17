%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f = - k^2 V exp(i k x)
%
% Solved directly using skeletonization with proxy surfaces
%
%%%%%

L = 10; % length of grid
N1 = 400; % number of grid points

zk = 10;

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

[coefs,dinds] = bump(xxgrid,yygrid,1.25,0.4,1,1e-12);
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

% Constructing identity + sparse corrections

[inds,corrs] = get_correct_helm(h,zk);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
idspmat = id_plus_corr_sum(coefs,spmats,dinds,h);

% Defining integral operators 

srcinfo = []; srcinfo.r = [xxgrid(dinds) yygrid(dinds)].'; srcinfo.wts = h^2*ones(length(dinds),1);
targinfo = []; targinfo.r = srcinfo.r; 

gfunc = @(s,t) helmgreen1(zk,s.r,t.r);
Afun = @(i,j) kern_matgen(i,j,srcinfo,targinfo,coefs,gfunc,idspmat);

% Solve with FLAM

quads = srcinfo.wts;
pxyf = @(x,slf,nbr,l,ctr)  pxyfun_helm(x,slf,nbr,l,ctr,quads,zk,V);

rs          = srcinfo.r;
occ         = 2048;
rank_or_tol = 1E-8;
opts        = [];

start = tic;
F = rskelf(Afun,rs,occ,rank_or_tol,pxyf,opts);
t2 = toc(start);
fprintf('%5.2e s : time to factorize inverse (skel) \n',t2)

start = tic;
sol = rskelf_sv(F,rhs_vec);
t3 = toc(start);
fprintf('%5.2e s : time to solve (skel) \n',t3)

mu = zeros(size(xxgrid));
mu(dinds) = sol;

% Plot with FFT

[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
kerns = kernmat(src,targ,@(s,t) helmgreen1(zk,s,t),h);
kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
% usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
usca = usca(:,:,1);

utot = usca + uinc;
utot = reshape(utot,size(xxgrid));

figure(3); clf
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
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,coefs,0.1,0.1,zk,dinds,'helm');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)

return
