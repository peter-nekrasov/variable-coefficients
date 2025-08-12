%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves
%
%
%%%%%

L = 1000; % length of grid
N1 = 100; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

[coefs,dinds,pcoefs] = bump2(xxgrid,yygrid,2,50,8,1e-12);
[iinds,jinds] = ind2sub(size(xxgrid),dinds);

a0 = pcoefs{1}; 
b0 = pcoefs{3}; 
g0 = pcoefs{5}; 

fprintf('Number of points: %d \n',size(dinds,1))

% Finding positive real roots
[rts,ejs] = find_roots(b0 / a0, g0 / a0);
ejs = ejs/a0;
zk = rts((abs(angle(rts)) < 1e-6) & (real(rts) > 0));

% RHS (Incident field)
theta = -pi/4;
pws = planewave(zk,[xxgrid(:) yygrid(:)].',theta,10);
phizinc = pws(:,:,1);
phiinc = pws(:,:,1) / zk;

uincs = zeros([size(dinds) 8]);
uincs(:,:,1) = pws(dinds,:,1);
uincs(:,:,2:4) = pws(dinds,:,4:6);
uincs(:,:,5) = pws(dinds,:,4) + pws(dinds,:,6);
uincs(:,:,6) = pws(dinds,:,7) + pws(dinds,:,9);
uincs(:,:,7) = pws(dinds,:,8) + pws(dinds,:,10);
uincs(:,:,8) = pws(dinds,:,1) / zk;

rhs_vec = get_rhs(coefs,uincs);
coefs(:,:,1:end-1) = 1/2*coefs(:,:,1:end-1);

figure(1); clf
tiledlayout(1,3);

nexttile
aplot = pcoefs{1} + pcoefs{2};
pcolor(xxgrid,yygrid,aplot); shading interp;
colorbar
title('\alpha')
drawnow

nexttile
bplot = pcoefs{3} + pcoefs{4};
pcolor(xxgrid,yygrid,bplot); shading interp;
colorbar
title('\beta')
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
[inds,corrs] = get_correct_fg(h,a0);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
idspmats = id_plus_corr_sum(coefs,spmats,dinds,h);

% Defining integral operators 

srcinfo = []; srcinfo.r = [xxgrid(dinds) yygrid(dinds)].'; srcinfo.wts = h^2*ones(length(dinds),1);
targinfo = []; targinfo.r = srcinfo.r; 

gfunc = @(s,t) fggreen(s.r,t.r,rts,ejs);
Afun = @(i,j) kern_matgen(i,j,srcinfo,targinfo,coefs,gfunc,idspmats);

% Solve with FLAM

quads = srcinfo.wts;
pxyf = @(x,slf,nbr,l,ctr)  pxyfun_fg(x,slf,nbr,l,ctr,quads,zk,rts,ejs,coefs);

rs          = srcinfo.r;
occ         = 2048;
rank_or_tol = 1E-8;
opts        = [];
% pxyf        = [];

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
gfunc = @(s,t) fggreen(s,t,rts,ejs);
kerns = kernmat(src,targ,gfunc,h);
kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = kerns(:,:,[1 8]);
evalspmats = {spmats{1},spmats{8}};

usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
% usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
phizsca = usca(:,:,1)/2;
phisca = usca(:,:,2);

phitot = phisca + phiinc;
phiztot = phizsca + phizinc;

phitot = reshape(phitot,size(xxgrid));
phiztot = reshape(phiztot,size(xxgrid));

figure(2);
pc = pcolor(xxgrid,yygrid,real(mu)); shading interp;
title('Re(\mu)')
colorbar

figure(3); clf
tiledlayout(2,2)

nexttile
pc = pcolor(xxgrid,yygrid,real(phitot)); shading interp;
title('Re(\phi)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(phitot)); shading interp;
title('|\phi|')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(phiztot)); shading interp;
title('real(\phi_z)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(phiztot)); shading interp;
title('|\phi_z|')
colorbar
       
% Calculate error with finite difference
utots = cat(3,phiztot,phitot);
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utots,h,pcoefs,10,10,zk,dinds,'fg');

fprintf('Absolute error (fin diff): %.4e \n',abs_err)
fprintf('Relative error (fin diff): %.4e \n',rel_err)

return