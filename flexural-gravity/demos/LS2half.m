%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural-gravity waves
%
%
%%%%%

L = 10; % length of grid
N1 = 300; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

a0 = 1.2; 
b0 = a0^(1/3); 
g0 = -4.4; 
nu = 0.2;

[acoef,dinds] = bump(xxgrid,yygrid,0.5,0.5,6,1e-12);
[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

ac = acoef(:,:,1);
ax = acoef(:,:,2);
ay = acoef(:,:,3);
axx = acoef(:,:,4);
axy = acoef(:,:,5);
ayy = acoef(:,:,6);
alpha = a0 + ac;
beta = (alpha).^(1/3);
bc = beta - b0;

coefs = zeros([size(dinds) 8]);
coefs(:,:,1) = (ac*b0 - a0*bc) / a0;
coefs(:,:,2) = -(1-nu)*ayy;
coefs(:,:,3) = 2*(1-nu)*axy;
coefs(:,:,4) = -(1-nu)*axx;
coefs(:,:,5) = (axx+ayy);
coefs(:,:,6) = 2*ax;
coefs(:,:,7) = 2*ay;
coefs(:,:,8) = -g0*ac/a0;
coefs = coefs*a0 ./ alpha;

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

figure(1); clf
tiledlayout(1,3);

nexttile
aplot = a0 + zeros(length(xxgrid(:)),1);
aplot(dinds) = alpha;
aplot = reshape(aplot,size(xxgrid));
pcolor(xxgrid,yygrid,aplot); shading interp;
colorbar
title('\alpha')
drawnow

nexttile
bplot = b0 + zeros(length(xxgrid(:)),1);
bplot(dinds) = beta;
bplot = reshape(bplot,size(xxgrid));
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
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
[inds,corrs] = get_correct_fg(h,a0);
gfunc = @(s,t) fggreen(s,t,rts,ejs);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
% kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);

% Solve with GMRES
start = tic;
sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,spmats,h,dinds,iinds,jinds,N2),rhs_vec,[],1e-12,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

evalkerns = zeros([N2 N2 2]);
evalkerns(:,:,1) = kerns(:,:,1);
evalkerns(:,:,2) = kerns(:,:,8);
evalspmats = {spmats{1},spmats{8}};

usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
phisca = usca(:,:,1);
phizsca = usca(:,:,2);

phitot = phisca + phiinc;
phiztot = phizsca + phizinc;

phitot = reshape(phitot,size(xxgrid));
phiztot = reshape(phiztot,size(xxgrid));

%%

figure(2);
tiledlayout(2,3)

nexttile
pc = pcolor(xxgrid,yygrid,real(mu)); shading interp;
title('Re(\mu)')
colorbar

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
title('real(\phi_n)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(phiztot)); shading interp;
title('|\phi_n|')
colorbar
       
% Calculate error with finite difference
utots = cat(3,phiztot,phitot);
err = get_fin_diff_err(xxgrid,yygrid,utots,h,coefs,0.1,0.1,zk,dinds,'fg');

fprintf('Finite difference error: %.4e \n',err)

return

%%

figure(1);
s = surf(X,Y,H);
s.EdgeColor = 'none';

figure(2);
s = surf(X,Y,real(phi_n));
s.EdgeColor = 'none';
