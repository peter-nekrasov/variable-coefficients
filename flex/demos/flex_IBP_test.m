% solve the forward problem given G and corrections

%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural waves
%
%%%%%

L = 1000; % length of grid
N1 = 101; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = ndgrid(xs);
h = xs(2) - xs(1);

 w = 8;
[~,~,pcoefs] = bump4(xxgrid,yygrid,0.5,50,w,1e-12);

a0 = pcoefs{1}; 
b0 = pcoefs{3}; 
zk = (b0/a0)^0.25;
clear pcoefs

% Constructing integral operators
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);
[inds,corrs] = get_correct_flex(h);
gfunc = @(s,t) flexgreen2(zk,s,t);
%kerns = kernmat(src,targ,gfunc,h);
kerns = kernmat(src,targ,gfunc,h,inds,corrs);
kerns = gen_fft_kerns(kerns,sz,ind);
kerns = kerns / a0;

amp = 2; width = 35;
amp = 2;
[coefs0,dinds0,pcoefs0,~,alpha0] = bump4(xxgrid,yygrid,amp,width,w,1e-12);
[iinds,jinds] = ind2sub(size(xxgrid),dinds0);


nsrc = 1;
f = (1:nsrc);
thetas = 2*pi*(1:nsrc)/nsrc;
srcs = L/4*[cos(thetas);sin(thetas)];
targ2 = [xxgrid(:) yygrid(:)].';
sensfuns = kernmat(srcs,targ2,gfunc,h);
% sensfuns = gfunc(srcs,targ2)*h*h;

sensfuns = f*sensfuns(:,:,1).'/a0;


fprintf('Number of points: %d \n',size(dinds0,1))

% RHS (Incident field)

% % incfield = 'pointsource';

src = 1*[-L/4;-L/4];
src = srcs;
uincs = flexgreen2(zk,src,[xxgrid(:) yygrid(:)].');

uinc = uincs(:,:,1);

rhs_vec0 = get_rhs(coefs0,uincs(dinds0,:,:));


coefspad = zeros(numel(xxgrid),1,size(coefs0,3));
coefspad(dinds0,:,:) = coefs0;

% vas = zeros(1,2);
% vals(1) = sum(alpha0(dinds0).*rhs_vec0.*conj(freshu_adj(dinds0,1))) * h^2;
% 
% rhs2 = get_rhs(coefs0,freshu_adj(dinds0,:,:))-freshadj_sol;
% freshu_adj2 = -alpha0(dinds0).*get_rhs(coefs0,freshu_adj(dinds0,:,:));
% 
% % vals(2) = sum(utot(dinds0).*conj(freshu_adj)) * h^2;
% vals(2) = sum(utot(dinds0).*conj(freshu_adj2)) * h^2;

width2 = 40;
vals =0;
% shift = [0;-L/4];
shift = [-L/8;0];
f = exp(-vecnorm(targ2(:,dinds0)+shift).^2/2/width2^2).'./alpha0(dinds0);
[usca0, sol0] = solve_flex(f, coefs0, kerns, iinds, jinds, N1, N2);

shift = [0;0.2*L];
g = exp(-vecnorm(targ2(:,dinds0)+shift).^2/2/width2^2).'./alpha0(dinds0);
[usca1, sol1] = solve_flex(g, coefs0, (kerns), iinds, jinds, N1, N2);

k = 1;
vals(1) = sum(alpha0(dinds0).^k.*f.*(usca1(dinds0,1))) * h^2;
vals(2) = sum(alpha0(dinds0).^k.*usca0(dinds0,1).*(g)) * h^2;

% vals(1) = sum(f.*conj(usca0(dinds0,1))) * h^2;
% vals(2) = sum(usca1(dinds0).*conj(g)) * h^2;
vals
abs(diff(vals)./vals(1:end-1))

% usens2 = sensfuns(:,dinds)*sol;
% freshsens = sensfuns(:,dinds0)*freshsol;

% [abs(usens2 - (usens)), abs(usens2 - (usens+delta*freshsens))]/abs(usens)

%%
figure(4);clf
subplot(1,2,1)
a = zeros(size(xxgrid));
a(dinds0) = log10(f);
% a = reshape(usca0(:,:,1),size(xxgrid));
pc = pcolor(xxgrid,yygrid,real(a)); shading interp;
clim([-15,inf])
colorbar
subplot(1,2,2)
a = zeros(size(xxgrid));
a(dinds0) = log10(g);
% a = reshape(usca1(:,:,1),size(xxgrid));
pc = pcolor(xxgrid,yygrid,real(a)); shading interp;
colorbar
clim([-15,inf])

return


%%


n = size(usca0,1);

tic;
Amat = zeros(n,n);
for i = 1:n

    x = zeros(n,1);
    x(i) = 1;
    src = targ2(:,i);
    % u = kernmat(src,targ2,gfunc,h,inds,corrs);
    u = gfunc(src,targ2);
    rhs_vec0 = get_rhs(coefs0,u(dinds0,:,:));
    rhs_vec0(i) = rhs_vec0(i) + 1;
    Amat(:,i) = rhs_vec0;
end

toc;


%%
i = floor((N1/2)^2 + N1/2);
figure(5);clf
subplot(1,2,1)
a = zeros(size(xxgrid));
a = reshape(Amat(:,i),size(xxgrid));
pc = pcolor(xxgrid,yygrid,real(a)); shading interp;

colorbar
subplot(1,2,2)
a = zeros(size(xxgrid));
a = reshape(Amat(i,:),size(xxgrid));
pc = pcolor(xxgrid,yygrid,real(a)); shading interp;
colorbar


function [usca, sol] = solve_flex(rhs_vec, coefs, kerns, iinds, jinds, N1, N2)
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,200);
evalkerns = kerns(:,:,1);
evalkerns = kerns;
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
end