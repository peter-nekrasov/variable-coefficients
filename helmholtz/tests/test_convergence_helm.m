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

L = 10;
Ns = 50:50:700; % needs to be an odd number for FFT

zk = 5;

errs = zeros(length(Ns),1);

for ii = 1:length(Ns)

N1 = Ns(ii);
xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);

h = xs(2) - xs(1);

[coefs,dinds] = bump(xxgrid,yygrid,-0.5,0.5,1,1e-12);
V = coefs(:,:,1);
coefs = coefs*zk^2;

[iinds,jinds] = ind2sub(size(xxgrid),dinds);

% RHS (Incident field)
theta = -pi/3;
uinc = planewave(zk,[xxgrid(:) yygrid(:)].',theta);
rhs_vec = get_rhs(coefs,uinc,dinds);

% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,0);
[inds,corrs] = get_correct_helm(h,zk);
gfunc = @(s,t) helmgreen1(zk,s,t);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
% kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);

% Solve with GMRES
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,spmats,h,dinds,iinds,jinds,N2),rhs_vec,[],1e-10,200);
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

% usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
usca = usca(:,:,1);

utot = usca + uinc;
utot = reshape(utot,size(xxgrid));

errs(ii) = get_fin_diff_err(xxgrid,yygrid,utot,h,coefs,0.1,0.1,zk,dinds,'helm');

end

%%

figure(2);
plot(log10(Ns),log10(errs),'x-');
hold on
plot(log10(Ns),log10(1e11*Ns.^(-8)))
xlabel('log_{10}(N)')
legend('error','N^{-8}')