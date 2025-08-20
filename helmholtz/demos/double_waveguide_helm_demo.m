%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f corresponds to the incident field of a Gaussian beam
%
% Solved directly using skeletonization
%
%%%%%

% Domain parameters

zk = 10;

L = 10;
N1 = 400; 

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

% Waveguide parameters

wgdist = 2+(2*pi/zk);
wglen = 0.9*L;
wgamp = (3.5^2-1);
wgwid = 0.1*L;

[coefs,dinds] = double_waveguide(xxgrid,yygrid,L,wgamp,wglen,wgwid,wgdist,0.08,1e-8);
V = coefs(:,:,1);
coefs = coefs*zk^2;

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

srcinfo = []; srcinfo.r = [xxgrid(dinds) yygrid(dinds)].'; srcinfo.wts = h^2*ones(length(dinds),1);
targinfo = []; targinfo.r = srcinfo.r; 

% RHS (Incident field)
k1 = zk;
k2 = 0;
src = []; src.r = [-L/2;-1/2*wgdist] + 1.5i*[1;0];
targ = []; targ.r = [xxgrid(:) yygrid(:)].';
uinc = helmgreen1(zk,src.r,targ.r);
uinc = uinc / max(abs(uinc(:)));
rhs_vec = get_rhs(coefs,uinc,dinds);

% RHS (Incident field)

figure(1); clf
tiledlayout(1,3);

nexttile
Vplot = zeros(length(xxgrid(:)),1);
Vplot(dinds) = V;
Vplot = reshape(Vplot,size(xxgrid));
pcolor(xxgrid,yygrid,Vplot); shading interp;
colorbar
title('V')
axis square
drawnow

nexttile
uincplot = reshape(uinc,size(xxgrid));
pcolor(xxgrid,yygrid,real(uincplot)); shading interp;
colorbar
title('u^{inc}')
axis square
drawnow

nexttile
rhsplot = zeros(length(xxgrid(:)),1);
rhsplot(dinds) = rhs_vec;
rhsplot = reshape(rhsplot,size(xxgrid));
pcolor(xxgrid,yygrid,real(rhsplot)); shading interp;
colorbar
title('rhs')
axis square
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

%% Plot with FFT

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

figure(2); clf
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
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,coefs,L/4,wgdist/2,zk,dinds,'helm');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)

return

%%% solving the waveguide example with FFT + GMRES (takes way longer)

gfunc = @(s,t) helmgreen1(zk,s,t);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);

% Solve with GMRES
start = tic;
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,spmats,h,dinds,iinds,jinds,N2),rhs_vec,[],1e-10,200);
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-8,5000);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('Time to solve: %5.2e s\n',t1)
