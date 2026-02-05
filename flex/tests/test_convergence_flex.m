L = 1000;
Ns = [51,101,201,301,401,501,601];

errs = zeros(length(Ns),1);
mus = zeros(length(Ns),1);

for ii = 1:length(Ns)

N1 = Ns(ii);
xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

[coefs,dinds,pcoefs] = bump4(xxgrid,yygrid,0.5,50,4,1e-12);
[iinds,jinds] = ind2sub(size(xxgrid),dinds);

a0 = pcoefs{1}; 
b0 = pcoefs{3}; 

fprintf('Number of points: %d \n',size(dinds,1))

% Finding positive real roots
zk = (b0/a0)^0.25;

theta = pi/4;
pws = planewave(zk,[xxgrid(:) yygrid(:)].',theta,10);
uinc = pws(:,:,1);
uincs = zeros([size(dinds) 7]);
uincs(:,:,1) = pws(dinds,:,1);
uincs(:,:,2:4) = pws(dinds,:,4:6);
uincs(:,:,5) = pws(dinds,:,4) + pws(dinds,:,6);
uincs(:,:,6) = pws(dinds,:,7) + pws(dinds,:,9);
uincs(:,:,7) = pws(dinds,:,8) + pws(dinds,:,10);
% uincs(:,:,8) = pws(dinds,:,1) / zk;

rhs_vec = get_rhs(coefs,uincs);

% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
[inds,corrs] = get_correct_flex(h);
gfunc = @(s,t) flexgreen2(zk,s,t);
spmats = get_sparse_corr_mat([N1 N1],inds,corrs);
idspmats = id_plus_corr_sum(coefs,spmats,dinds,h);
%kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);
kerns = kerns / a0;

% Solve with GMRES
start = tic;
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,200);
% sol = gmres(@(mu) fast_apply_fft_sub(mu,kerns,coefs,idspmats,iinds,jinds,N2),rhs_vec,[],1e-10,200);
mu = zeros(size(xxgrid));
mu(dinds) = sol;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)

evalkerns = kerns(:,:,1);
evalspmats = {spmats{1}};

% usca = sol_eval_fft_sub(sol,evalkerns,evalspmats,h,dinds,iinds,jinds,N1,N2);
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);

utot = uinc + usca;

utot = reshape(utot,size(xxgrid));
       
mus(ii) = mu((N1-1)/2+1,(N1-1)/2+1);
errs(ii) =  get_fin_diff_err(xxgrid,yygrid,utot,h,pcoefs,20,20,zk,dinds,'flex');

end
mus = abs(mus(1:end-1) - mus(end)) / abs(mus(end));


%%


figure(2); clf
t = tiledlayout(1,2);
nexttile
plot(log10(Ns(1:end-1)),log10(mus),'x-');
hold on
plot(log10(Ns),log10(1e6*Ns.^(-6)))
xlabel('log_{10}(N)')
legend('error','N^{-6}')
title('self-convergence')

nexttile
plot(log10(Ns),log10(errs),'x-');
hold on
plot(log10(Ns),log10(1e3*Ns.^(-6)))
xlabel('log_{10}(N)')
legend('error','N^{-6}')
title('finite difference')