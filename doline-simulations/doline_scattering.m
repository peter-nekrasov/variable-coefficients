ws = linspace(0.01,20,1000);

for ii = 1:1000
    doline_scattering_1(ws(ii),ii);
end

function doline_scattering_1(w,ii)

L = 2000; % length of grid
N1 = 400; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

[coefs,dinds,pcoefs,H] = doline(xxgrid,yygrid,-200,130,330,w,1e-12);
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

incfield = 'pointsource';
% incfield = 'planewave';

if strcmpi(incfield,'planewave')

theta = pi/4;
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

elseif strcmpi(incfield,'pointsource')

uincs = fggreen([-L/2;-L/2],[xxgrid(:) yygrid(:)].',rts,ejs);
phizinc = uincs(:,:,1)/2;
phiinc = uincs(:,:,end);
uincs = uincs(dinds,:,:);
uincs(:,:,1:end-1) = uincs(:,:,1:end-1)/2;

end

rhs_vec = get_rhs(coefs,uincs);
coefs(:,:,1:end-1) = coefs(:,:,1:end-1)/2;

% figure(1); clf
% tiledlayout(2,2);
% 
% nexttile
% aplot = pcoefs{1} + pcoefs{2};
% pcolor(xxgrid,yygrid,aplot); shading interp;
% colorbar
% title('\alpha')
% drawnow
% 
% nexttile
% bplot = pcoefs{3} + pcoefs{4};
% pcolor(xxgrid,yygrid,bplot); shading interp;
% colorbar
% title('\beta')
% drawnow
% 
% nexttile
% phizplot = reshape(phizinc,size(xxgrid));
% pcolor(xxgrid,yygrid,real(phizplot)); shading interp;
% colorbar
% title('\phi^{(i)}')
% drawnow
% 
% nexttile
% rhsplot = zeros(length(xxgrid(:)),1);
% rhsplot(dinds) = rhs_vec;
% rhsplot = reshape(rhsplot,size(xxgrid));
% pcolor(xxgrid,yygrid,real(rhsplot)); shading interp;
% colorbar
% title('rhs')
% drawnow
% 
% figure(2); clf
% pcolor(xxgrid,yygrid,H,'EdgeColor','none')
% title("Thickness (H)")
% colorbar 


% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
[inds,corrs] = get_correct_fg(h,a0);
gfunc = @(s,t) fggreen(s,t,rts,ejs);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
idspmats = id_plus_corr_sum(coefs,spmats,dinds,h);
%kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);

% Solve with GMRES
start = tic;
[sol,flag,relres,iter,resvec] = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-8,2000);
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,idspmats,iinds,jinds,N2),rhs_vec,[],1e-10,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

evalkerns = kerns(:,:,[1,8]);
evalspmats = {spmats{1},spmats{8}};

% usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
phizsca = usca(:,:,1)/2;
phisca = usca(:,:,2);

phiztot = phizsca + phizinc;
phitot = phisca + phiinc;

phiztot = reshape(phiztot,size(xxgrid));
phitot = reshape(phitot,size(xxgrid));

% figure(2);
% pc = pcolor(xxgrid,yygrid,real(mu)); shading interp;
% title('Re(\mu)')
% colorbar
% 
% figure(3); clf
% tiledlayout(2,2)
% 
% nexttile
% pc = pcolor(xxgrid,yygrid,real(phitot)); shading interp;
% title('Re(\phi)')
% colorbar
% 
% nexttile
% pc = pcolor(xxgrid,yygrid,abs(phitot)); shading interp;
% title('|\phi|')
% colorbar
% 
% nexttile
% pc = pcolor(xxgrid,yygrid,real(phiztot)); shading interp;
% title('real(\phi_z)')
% colorbar
% 
% nexttile
% pc = pcolor(xxgrid,yygrid,abs(phiztot)); shading interp;
% title('|\phi_z|')
% colorbar
       
% Calculate error with finite difference
utots = cat(3,phiztot,phitot);
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utots,h,pcoefs,10,10,zk,dinds,'fg');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)

save("doline_scattering_" + ii + ".mat", "sol","mu","w","ii","phizsca","phisca","flag","relres","iter","resvec","a0","b0","g0","pcoefs","xxgrid","yygrid","H","phiztot","phitot")

end