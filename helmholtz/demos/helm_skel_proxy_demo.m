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

L = 10;
N = 201; 

zk = 8;

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

[coefs,dinds] = bump(xxgrid,yygrid,-0.5,0.5,1,1e-12);
V = coefs{1};
coefs{1} = coefs{1}*zk^2;

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

% RHS (Incident field)
theta = -pi/3;
uinc = planewave(zk,[xxgrid(:) yygrid(:)].',theta);
uincs = {uinc};
rhs_vec = get_rhs(coefs,uincs,dinds);

figure(1); clf
tiledlayout(1,2);

nexttile
Vplot = zeros(length(xxgrid(:)),1);
Vplot(dinds) = V;
Vplot = reshape(Vplot,size(xxgrid));
s = pcolor(xxgrid,yygrid,Vplot);
s.EdgeColor = 'None';
colorbar
title('V')
drawnow

nexttile
rhsplot = zeros(length(xxgrid(:)),1);
rhsplot(dinds) = rhs_vec;
rhsplot = reshape(rhsplot,size(xxgrid));
s = pcolor(xxgrid,yygrid,real(rhsplot));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow

% Constructing identity + sparse corrections

[inds,corrs] = get_correct_helm(h,zk);
spmats = get_sparse_corr_mat(size(xxgrid),inds,corrs);
idspmat = id_plus_corr_sum(coefs,spmats,dinds,h);

% Defining integral operators 

gfunc = @(s,t) helmgreen1(zk,s,t);
Afun = @(i,j) kern_matgen(i,j,srcinfo,targinfo,idspmat,kernfun);

% Solve with FLAM

quads = srcinfo.wts;

pxyf = @(x,slf,nbr,l,ctr)  pxyfun_helm(x,slf,nbr,l,ctr,quads,zk,targinfo.V);

rs          = srcinfo.r;
occ         = 4000;
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

[src,targ,ind,sz] = get_fft_grid(N,L);
kerns = kernmat(src,targ,@(s,t) helm2d.green_cell_helm(zk,s,t),h);
kerns = gen_fft_kerns2(kerns,sz,ind);

evalkerns = {kerns{1}};
evalcorrs = {spmats{1}};

usca = sol_eval_fft_sub_helm(sol,evalkerns,evalcorrs,h,dinds,iinds,jinds,xxgrid);

utot = usca + uinc;

figure(2); clf
tiledlayout(1,3)

nexttile
pc = pcolor(xxgrid,yygrid,real(mu));
pc.EdgeColor = 'none';
title('Re(\mu)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(utot));
pc.EdgeColor = 'none';
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot));
pc.EdgeColor = 'none';
title('|\phi|')
colorbar
       
% Calculate error with finite difference
err = get_fin_diff_err_helm(xxgrid,yygrid,utot,h,coefs,0.1,0.1,zk);

fprintf('Finite difference error: %.4e \n',err)

return
