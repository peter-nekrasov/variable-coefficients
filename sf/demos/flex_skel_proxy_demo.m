%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% Helmholtz scattering problem:
%
%      \Delta^2 u - k^4 (1 + V(x) ) u = f
%
% in this case, f either corresponds to a plane wave or a point source
%
% Solved iteratively using FFT + GMRES 
%
%%%%%

% incfield = 'pointsource';
incfield = 'planewave';

L = 20; % length of grid
N1 = 201; % number of grid points

zk = 2;

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

[coefs,dinds] = bump(xxgrid,yygrid,1.5,1,1,1e-8);
V = coefs(:,:,1);
coefs = -coefs*zk^4;

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

% RHS (Incident field)

if strcmpi(incfield,'planewave')

theta = -pi/4;
uinc = planewave(zk,[xxgrid(:) yygrid(:)].',theta);

elseif strcmpi(incfield,'pointsource')

uinc = flexgreen(zk,[-L/3;-L/3],[xxgrid(:) yygrid(:)].');

end

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
[inds,corrs] = get_correct_flex_simplified(h);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
idspmat = id_plus_corr_sum(coefs,spmats,dinds,h);

% Defining integral operators

srcinfo = []; srcinfo.r = [xxgrid(dinds) yygrid(dinds)].'; srcinfo.wts = h^2*ones(length(dinds),1);
targinfo = []; targinfo.r = srcinfo.r; 

gfunc = @(s,t) flexgreen(zk,s.r,t.r);
Afun = @(i,j) kern_matgen(i,j,srcinfo,targinfo,coefs,gfunc,idspmat);

% Solve with FLAM

quads = srcinfo.wts;
pxyf = @(x,slf,nbr,l,ctr)  pxyfun_flex_simplified(x,slf,nbr,l,ctr,quads,zk,V);

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
kerns = kernmat(src,targ,@(s,t) flexgreen(zk,s,t),h);
kerns = gen_fft_kerns(kerns,sz,ind);

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
% usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
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
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,coefs,0.1,0.1,zk,dinds,'flex-simplified');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)

return
