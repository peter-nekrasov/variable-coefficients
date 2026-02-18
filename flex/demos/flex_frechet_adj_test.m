% solve the forward problem given G and corrections

%%%%%
%
% Solving the adjointed Lippman-Schwinger equation for plane wave 
% scattering of flexural waves
%
%%%%%

L = 1000; % length of grid
N1 = 201; % number of grid points

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
amp = 0.5;
[coefs0,dinds0,pcoefs0,~,alpha0] = bump4(xxgrid,yygrid,amp,width,w,1e-12);
[iinds,jinds] = ind2sub(size(xxgrid),dinds0);


nsrc = 1;
f = (1:nsrc);
thetas = 2*pi*(1:nsrc)/nsrc;
srcs = L*[cos(thetas);sin(thetas)];
targ2 = [xxgrid(:) yygrid(:)].';
sensfuns = kernmat(srcs,targ2,gfunc,h);
% sensfuns = gfunc(srcs,targ2)*h*h;

sensfuns = f*sensfuns(:,:,1).'/a0;


fprintf('Number of points: %d \n',size(dinds0,1))

% RHS (Incident field)

% % incfield = 'pointsource';

% uincs = flexgreen2(zk,[-L/2;-L/2],[xxgrid(:) yygrid(:)].');

theta =0;
pws = planewave(zk,[xxgrid(:) yygrid(:)].',theta,10);

uincs = zeros([size(xxgrid(:)) 7]);
uincs(:,:,1) = pws(:,:,1);
uincs(:,:,2:4) = pws(:,:,4:6);
uincs(:,:,5) = pws(:,:,4) + pws(:,:,6);
uincs(:,:,6) = pws(:,:,7) + pws(:,:,9);
uincs(:,:,7) = pws(:,:,8) + pws(:,:,10);

uinc = uincs(:,:,1);

rhs_vec0 = get_rhs(coefs0,uincs(dinds0,:,:));


start = tic;
[usca0, sol0] = solve_flex(rhs_vec0, coefs0, kerns, iinds, jinds, N1, N2);
mu = zeros(size(xxgrid));
mu(dinds0) = sol0;
t1 = toc(start);
fprintf('%5.2e s : time to solve\n',t1)


usens = sensfuns(:,dinds0)*sol0;

coefspad = zeros(numel(xxgrid),1,size(coefs0,3));
coefspad(dinds0,:,:) = coefs0;

utot = uinc + usca0(:,:,1);
utots = uincs + usca0;

deltas = 10.^(-3:-1:-4);
dus = zeros(2,length(deltas));
vals = zeros(2,length(deltas));
for i = 1:length(deltas)
delta = deltas(i);
[coefs,dinds,~,~,alpha] = bump4(xxgrid,yygrid,amp+delta,width,w,1e-12);
rhs_vec = get_rhs(coefs,uincs(dinds,:,:));
[iinds,jinds] = ind2sub(size(xxgrid),dinds);
[usca, sol] = solve_flex(rhs_vec, coefs, kerns, iinds, jinds, N1, N2);


[coefsp,dindsp,pcoefsp,~,alphap] = bump4(xxgrid,yygrid,amp+delta,width,w,1e-12);
[coefsm,dindsm,pcoefsm,~,alpham] = bump4(xxgrid,yygrid,amp-delta,width,w,1e-12);
coefsmpad = zeros(numel(alpham(:)),1,size(coefs,3));
coefsmpad(dindsm,:,:) = -coefsm.*alpham(dindsm);
coefsmpad(dindsp,:,:) = coefsmpad(dindsp,:,:)+ coefsp.*alphap(dindsp);
coefsmpad = coefsmpad/2/delta;
dcoefs = coefsmpad(dinds0,:,:);
dalpha = (alphap-alpham)/2/delta./alpha0;
dcoefs = dcoefs./alpha0(dinds0);

[iinds,jinds] = ind2sub(size(xxgrid),dinds0);

% dcoefs(:,:,[1:4]) = 0;
rhsfresh = get_rhs(dcoefs,utots(dinds0,:,:)) - (dalpha(dinds0)).*mu(dinds0);

% rhsfresh = fast_apply_fft(mu(dinds0),kerns,dcoefs,iinds,jinds,N2);
% rhsfresh = rhsfresh - mu(dinds0);
% rhsfresh = -(rhsfresh + (dalpha(dinds0)).*mu(dinds0));
% rhsfresh = rhsfresh + get_rhs(dcoefs,uincs(dinds0,:,:));

[freshu, freshsol] = solve_flex(rhsfresh, coefs0, kerns, iinds, jinds, N1, N2);
dus(1,i) = norm(usca(:,:,1) - usca0(:,:,1)) / norm(usca(:,:,1));
dus(2,i) = norm(usca(:,:,1) - freshu(:,:,1)*delta - usca0(:,:,1)) / norm(usca(:,:,1));

vals(1,i) = sensfuns(:,dinds0)*freshsol;

usrcs = flexgreen2(zk,srcs,[xxgrid(:) yygrid(:)].');
usrcs = conj(sum(usrcs.*f,2));

rhs_adj = (get_rhs(coefs0,usrcs(dinds0,:,:)));
[freshu_adj, freshadj_sol] = solve_flex(rhs_adj, coefs0, conj(kerns), iinds, jinds, N1, N2);
freshu_adj = freshu_adj + (usrcs(:,:,:));

vals(2,i) = sum(alpha0(dinds0).*rhsfresh.*conj(freshu_adj(dinds0,1))) * h^2;

freshu_adj = alpha0(dinds0).*(get_rhs(dcoefs,freshu_adj(dinds0,:,:)) - (dalpha(dinds0)).*freshadj_sol);

vals(3,i) = sum(utot(dinds0).*conj(freshu_adj)) * h^2;




end


[deltas; dus]
vals
abs(diff(vals))

usens2 = sensfuns(:,dinds)*sol;
freshsens = sensfuns(:,dinds0)*freshsol;

% [abs(usens2 - (usens)), abs(usens2 - (usens+delta*freshsens))]/abs(usens)

utot = reshape(utot,size(xxgrid));

figure(2);
pc = pcolor(xxgrid,yygrid,real(mu)); shading interp;
title('Re(\mu)')
colorbar

figure(3); clf
tiledlayout(1,2)

nexttile
pc = pcolor(xxgrid,yygrid,real(utot)); shading interp;
title('real(u)')
colorbar

nexttile
pc = pcolor(xxgrid,yygrid,abs(utot)); shading interp;
title('|u|')
colorbar
       


% %%
% figure(3);clf
% tiledlayout(1,2)
% nexttile
% pc = pcolor(xxgrid,yygrid,reshape(real(usca - usca0)/delta,size(xxgrid))); shading interp;
% colorbar
% % xlim([-10,10]*5)
% % ylim([-10,10]*5)       
% 
% nexttile
% upad = zeros(size(xxgrid)); upad(dinds) = rhsfresh;
% pc = pcolor(xxgrid,yygrid,reshape(real(freshu),size(xxgrid))); shading interp;
% % pc = pcolor(xxgrid,yygrid,abs(upad)); shading interp;
% colorbar
% % xlim([-10,10]*5)       
% % ylim([-10,10]*5)       

% Calculate error with finite difference
[abs_err,rel_err] = get_fin_diff_err(xxgrid,yygrid,utot,h,pcoefs0,10,10,zk,dinds0,'flex');

fprintf('Absolute error: %.4e \n',abs_err)
fprintf('Relative error: %.4e \n',rel_err)

return



function [usca, sol] = solve_flex(rhs_vec, coefs, kerns, iinds, jinds, N1, N2)
sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,200);
evalkerns = kerns(:,:,1);
evalkerns = kerns;
usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
end