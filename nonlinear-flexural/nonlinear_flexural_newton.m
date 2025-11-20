%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for the 
% nonlinear simplified flexural scattering problem:
%
%      \Delta^2 u - k^4 (1 + V(x) (a0 + a1 | u |^2 ) ) u = 0
%
% Solved iteratively using poor man's Newton
%
%%%%%

% incfield = 'pointsource';
% incfield = 'planewave';
incfield = 'beam';

a0 = 1;
a1 = 2;
zk = 4;

step_size = 0.1;

L = 10; % length of grid
N1 = 500; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

wgdist = 2+(2*pi/zk);
wglen = 0.8*L;
wgamp = (2^4-1); % 0.25
wgwid = 0.1*L;

[coefs,dinds] = double_waveguide(xxgrid,yygrid,L,wgamp,wglen,wgwid,wgdist,0.08,1e-8);
V = coefs(:,:,1);

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

fprintf('Number of points: %d \n',size(dinds,1))

% RHS (Incident field)

if strcmpi(incfield,'planewave')

theta = -pi/4;
uinc = planewave(zk,[xxgrid(:) yygrid(:)].',theta);

elseif strcmpi(incfield,'pointsource')

uinc = flexgreen(zk,[-L/3;-L/3],[xxgrid(:) yygrid(:)].');

elseif strcmpi(incfield,'beam')

k1 = zk;
k2 = 0;
src = []; src.r = [-L/2;-1/2*wgdist] + 4i*[1;0];
targ = []; targ.r = [xxgrid(:) yygrid(:)].';
uinc = helmgreen1(zk,src.r,targ.r) + helmgreen1(1i*zk,src.r,targ.r);
uinc = uinc / max(abs(uinc(:)));

end

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
uincplot = reshape(uinc,size(xxgrid));
pcolor(xxgrid,yygrid,real(uincplot)); shading interp;
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
utot_k = uinc(dinds,:,:);

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

% rhs_pde = zk.^4*coefs.*(a0 + a1*uinc(dinds,:,:).*conj(uinc(dinds,:,:))).*uinc(dinds,:,:);

tol = 1e-6;

a1_k = 0;

a1vec = [];
itervec = [];
rescell = {};

while a1_k < a1

fprintf("a1 = %f. \n",a1_k)

a1vec = [a1vec a1_k];

ie_resid = 1;
iter = 0;

resvec = mu_k - zk^4*V.*(a0 + a1_k*(utot_k).*conj(utot_k)).*utot_k;
resvec = vecnorm(resvec) / vecnorm(mu_k) ;

while (ie_resid) > tol

iter = iter + 1;

rhs_vec = - mu_k + zk^4*V.*(a0 + a1_k*(utot_k).*conj(utot_k)).*utot_k;
rhs_vec = [rhs_vec; conj(rhs_vec)];   

% rhs_vec = zk^4*V.*(a0 + a1_k*(utot_k).*conj(utot_k)).*uinc(dinds,:,:)  ;
% rhsnorm = vecnorm(rhs_vec);

coef1 = -zk.^4*V.*(a0+2*a1_k.*utot_k.*conj(utot_k));
coef2 = -zk.^4*V.*(a1_k.*utot_k.^2);

coefs = cat(3,coef1,coef2);

% Solve with GMRES
start = tic;
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,idspmats,iinds,jinds,N2),rhs_vec,[],1e-10,200);
[delta_mu_k,flag,relres,itv,resv] = gmres(@(mu) fast_apply_fft_newton(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-4,5000);
t1 = toc(start);

delta_mu_k = delta_mu_k(1:numel(dinds));

mu_k = mu_k + delta_mu_k;

% usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
u_k = sol_eval_fft(mu_k,evalkerns,iinds,jinds,N1,N2);
u_k = u_k(dinds,:,1);
utot_k = u_k + uinc(dinds,:,:);

ie_resid = mu_k - zk^4*V.*(a0 + a1_k*(utot_k).*conj(utot_k)).*utot_k ;
ie_resid = vecnorm(ie_resid) / vecnorm(mu_k);

if flag == 0
    fprintf("Iteration %d. GMRES converged (%d). PDE Residual: %.3e \n",iter,itv(2),ie_resid)
else
    fprintf("Iteration %d. GMRES failed to converge (%d). PDE Residual: %.3e \n",iter,itv(2),ie_resid)
end

resvec = [resvec ie_resid];

if numel(resvec) > 2
    if resvec(end) > resvec(end-1) && resvec(end-1) > resvec(end-2)
        mu_k = mu_k_prev;
        utot_k = utot_k_prev;
        a1_k = a1_k_prev;
        step_size = step_size/2;
        break
    end
end


end

itervec = [itervec iter];
rescell{end+1} = resvec;

if a1_k == 0
    usca = sol_eval_fft(mu_k,evalkerns,iinds,jinds,N1,N2);
    mulin = mu_k;
    ulin = usca + uinc;
    ulin = reshape(ulin,size(xxgrid));
end

a1_k_prev = a1_k;
a1_k = a1_k + step_size;
mu_k_prev = mu_k;
utot_k_prev = utot_k;

end

%%

mu_tot = zeros(size(xxgrid));
mu_tot(dinds) = mu_k;

usca = sol_eval_fft(mu_k,evalkerns,iinds,jinds,N1,N2);

utot = usca + uinc;
utot = reshape(utot,size(xxgrid));

save("flex_nonlinear_waveguide_" + zk + "_" + a1 + ".mat","a0","a1","zk","a1_k","xxgrid","yygrid","coefs","V","uinc","a1vec","itervec","rescell","mu_k","usca","utot","mu_tot","ulin","mu_lin")


figure(2);
tiledlayout(2,2)

% nexttile
% pc = pcolor(xxgrid,yygrid,real(mu_tot)); shading interp;
% title('Re(\mu)')
% colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(ulin)); shading interp;
title('Re(u), \alpha_1 = 0')
clim([-0.9 0.9])
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(ulin)); shading interp;
title('|u|, \alpha_1 = 0')
clim([0 1.2])
colorbar

% nexttile
% pc = pcolor(xxgrid,yygrid,real(mu_tot)); shading interp;
% title('Re(\mu)')
% colorbar

nexttile
pc = pcolor(xxgrid,yygrid,real(utot)); shading interp;
title("Re(u), \alpha_1 = " + a1_k)
clim([-0.9 0.9])
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot)); shading interp;
title("|u|, \alpha_1 = " + a1_k)
clim([0 1.2])
colorbar
       
% Calculate error with finite difference

pdecoef = -zk^4*Vtot.*(a0 + a1_k*utot.*conj(utot));
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,pdecoef,1,-2,zk,dinds,'flex-nonlinear');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)