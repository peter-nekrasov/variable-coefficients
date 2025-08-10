%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta u + k^2 (1 + V(x) ) u = f
%
% in this case, f = - k^2 V exp(i k x)
%
% Solved directly using skeletonization
%
%%%%%

% Domain parameters

zk = 10;

L = 5;
N = 401; 

xs = L*(-floor(N/2):floor(N/2))/floor(N/2);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

% Waveguide parameters

wgdist = 2+(2*pi/zk);
wglen = 1.7*L;
wgamp = (3.5^2-1);
wgwid = 0.17*L;

coefs = double_waveguide(xxgrid,yygrid,wgamp,wglen,wgwid,wgdist,0.08);
V = coefs{1};

dinds = find(abs(V) > 1e-8 );
[iinds,jinds] = find(abs(V) > 1e-8 );

fprintf('Number of points: %d \n',size(dinds,1))

srcinfo = []; srcinfo.r = [xxgrid(dinds) yygrid(dinds)].'; srcinfo.wts = h^2*ones(length(dinds),1);
targinfo = []; targinfo.r = [xxgrid(dinds) yygrid(dinds)].'; 
targinfo.V = V(dinds);

% RHS (Incident field)
k1 = zk;
k2 = 0;
src = []; src.r = [-L;-1/2*wgdist] + 1.5i*[1;0];
targ = []; targ.r = [xxgrid(:) yygrid(:)].';
uinc = helm2d.green_cell_helm(zk,src.r,targ.r);
uinc = reshape(uinc{1},size(xxgrid)) / max(uinc{1}(:));
[rhs_vec, rhs] = get_rhs_vec_helm(coefs,zk,uinc);
rhs_vec = rhs_vec(dinds);

figure(1); clf
tiledlayout(1,3);

nexttile
s = pcolor(xxgrid,yygrid,V);
s.EdgeColor = 'None';
colorbar
title('V')
drawnow

nexttile
s = pcolor(xxgrid,yygrid,real(uinc));
s.EdgeColor = 'None';
colorbar
title('u^{inc}')
drawnow

nexttile
s = pcolor(xxgrid,yygrid,real(rhs));
s.EdgeColor = 'None';
colorbar
title('rhs')
drawnow
 

% Constructing identity + sparse corrections

[inds,corrs] = get_correct_helm(h,zk);
spmats = get_sparse_corr(size(xxgrid),inds,corrs);
idspmat = id_plus_corr_sum_helm(zk,coefs,spmats,dinds,h);

% Defining integral operators 

kernfun = @(s,t) kern_sum_helm(zk,s,t);
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
title('Re(u)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot));
pc.EdgeColor = 'none';
title('|u|')
colorbar
       

%%

figure(4); clf
theme(gcf,"light")

pc = pcolor(xxgrid,yygrid,abs(utot));
pc.EdgeColor = 'none';
title('|u|')
colorbar
axis square
       

% Calculate error with finite difference
err = get_fin_diff_err_helm(xxgrid,yygrid,utot,h,coefs,L/2,wgdist/2,zk);

fprintf('Finite difference error: %.4e \n',err)


return
