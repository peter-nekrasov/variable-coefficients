% Time domain simulation of a line of sources

% Spatial parameters

incfield = 'linesource';
srcloc = [-1500;0];
theta = 0;

N1 = 400;

x1 = -1E3;
x2 = 5E3;

y1 = -3E3;
y2 = 3E3;

if (x2 - x1) ~= (y2 - y1)
    error('horizontal and vertical lengths must be equal')
end

L = x2 - x1; 

% Frequency/time parameters

t0 = 0;
w0 = 0.725;
sigma = 0.12;

gtilde = @(w) exp(1i*t0*(w-w0)).*exp(-(w-w0).^2/(2*sigma.^2));
 
nw = 300;
freqs = linspace(w0 - 4*sigma,w0 + 4*sigma,nw);

xs = linspace(x1,x2,N1);
ys = linspace(y1,y2,N1);
[xxgrid,yygrid] = meshgrid(xs,ys);
h = xs(2) - xs(1);
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,0);

[coefs, dinds, pcoefs, H] = rolls(xxgrid,yygrid,0,40E2,-25E2,25E2,0.75,333.3,freqs(1),1e-8); % remove gbar from coefs vector
[iinds,jinds] = ind2sub(size(xxgrid),dinds);
a0 = pcoefs{1}; 
[inds,corrs] = get_correct_fg(h,a0);

mus = zeros(length(freqs), length(dinds));
phis = zeros(length(freqs), length(xxgrid(:)));
phizs = zeros(length(freqs), length(xxgrid(:)));
flags = zeros(length(freqs));
relreses = zeros(length(freqs));
iters = zeros(length(freqs),2);
resvecs = cell(length(freqs));

mu_pre = zeros(length(dinds),1);

for ii = 1:length(freqs)

    w = freqs(ii);
    fprintf('Frequency: %.4f ', w);

    [coefs,~, pcoefs] = rolls(xxgrid,yygrid,0,40E2,-25E2,25E2,0.75,333.3,w,1e-8); % remove gbar from coefs vector

    a0 = pcoefs{1};
    alpha = pcoefs{1} + pcoefs{2};
    b0 = pcoefs{3}; 
    g0 = pcoefs{5}; 

    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    zk = rts((imag(rts) == 0) & (real(rts) > 0));
    ejs = ejs/a0;

    if strcmpi(incfield,'linesource')

    uincs = green1d(rts,ejs,srcloc,[xxgrid(:) yygrid(:)].',theta);
    phizinc = uincs(:,:,1)/2;
    phiinc = uincs(:,:,end);
    uincs = uincs(dinds,:,:);
    uincs(:,:,1:end-1) = uincs(:,:,1:end-1)/2;
    
    elseif strcmpi(incfield,'pointsource')
    
    uincs = fggreen([x1;0],[xxgrid(:) yygrid(:)].',rts,ejs);
    phizinc = uincs(:,:,1)/2;
    phiinc = uincs(:,:,end);
    uincs = uincs(dinds,:,:);
    uincs(:,:,1:end-1) = uincs(:,:,1:end-1)/2;
    
    end

    rhs_vec = get_rhs(coefs,uincs);
    coefs(:,:,1:end-1) = 1/2*coefs(:,:,1:end-1);

    gfunc = @(s,t) fggreen(s,t,rts,ejs);
    kerns = kernmat(src,targ,gfunc,h,inds,corrs);
    kerns = gen_fft_kerns(kerns,sz,ind);
    
    evalkerns = kerns(:,:,[1 8]);

    % Solve with GMRES
    start = tic;
    % Solve with GMRES
    [mu, flag, relres, iter, resvec] = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2,alpha(dinds)/a0),rhs_vec,[],1e-6,5000,[],[],mu_pre);

    mu_pre = mu;
    if flag 
        fprintf('GMRES failed with relative residual %f. \n', relres);
    else 
        fprintf('GMRES succeeded in %d iterations. \n', iter(2));    
    end    
    
    sol = zeros(size(xxgrid));
    sol(dinds) = mu;
    t1 = toc(start);
    fprintf('%5.2e s : time to solve\n',t1)
    
    usca = sol_eval_fft(mu,evalkerns,iinds,jinds,N1,N2);
    phizsca = usca(:,:,1)/2;
    phisca = usca(:,:,2);
        
    phitot = phisca + phiinc;
    phiztot = phizsca + phizinc;

    mus(ii,:) = mu;
    phis(ii,:) = phitot(:);
    phizs(ii,:) = phiztot(:);
    flags(ii) = flag;
    relreses(ii) = relres;
    resvecs{ii} = resvec;

    f = figure(1);
    set(f,'Position',[1 1200 1332 300])
    t = tiledlayout(1,3);
    nexttile
    phiplot = reshape(phitot,size(xxgrid));
    pc = pcolor(xxgrid,yygrid,imag(phiplot));
    pc.EdgeColor = 'none';
    title('\Im(\phi)')
    colorbar

    nexttile
    phizplot = reshape(phiztot,size(xxgrid));
    pc = pcolor(xxgrid,yygrid,imag(phizplot));
    pc.EdgeColor = 'none';
    title('\Im(\partial_z\phi)')
    colorbar

    nexttile
    muplot = reshape(sol,size(xxgrid));
    pc = pcolor(xxgrid,yygrid,imag(muplot));
    pc.EdgeColor = 'none';
    title('Re(\mu)')
    colorbar
    title(t,['\omega = ',num2str(w)])
    drawnow

    utots = cat(3,reshape(phiztot,size(xxgrid)),reshape(phitot,size(xxgrid)));
    [abs_err, rel_err] = get_fin_diff_err(xxgrid,yygrid,utots,h,pcoefs,50,50,zk,dinds,'fg');

    fprintf('Absolute error (fin diff): %.4e \n',abs_err)
    fprintf('Relative error (fin diff): %.4e \n',rel_err)

    ints_phi(ii,:) = phitot.*gtilde(freqs(ii));
    ints_phi_z(ii,:) = phiztot.*gtilde(freqs(ii));
end

return

%%

figure(2); clf
plot(freqs,real(phizs(:,500)));
hold on
plot(freqs,imag(phizs(:,500)));
xlim([freqs(1) freqs(end)])

figure(3); clf 
plot(freqs,real(exp(1i*200*freqs).*phizs(:,500).'));
hold on
plot(freqs,imag(exp(1i*200*freqs).*phizs(:,500).'));
title('exp(200i\xi)u(x_0,\xi)')
xlabel('\xi')
xlim([freqs(1) freqs(end)])

%%

ntimes = 1000;
ts = linspace(-10,200,ntimes);
expmat = exp(-1i*ts(:)*freqs);

phi_sols = expmat*ints_phi;
phi_z_sols = expmat*ints_phi_z;

%%

v = VideoWriter('roll_scattering_4.mp4', 'MPEG-4');  
v.FrameRate = 25;  % adjust as desired
open(v);

f = figure('Position',[1 1 1288 544]); clf
f.Theme = 'light';
t = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');

nexttile(1)
zztarg = reshape(phi_z_sols(1,:),size(xxgrid));
h1 = pcolor(xxgrid, yygrid, imag(zztarg),'EdgeColor','interp');
shading interp;
hold on
contour(xxgrid,yygrid,H,[4 4],'k','ShowText','off',"LabelFormat","%0.1f m")
clim(0.8*[min(real(phi_z_sols(:))), max(real(phi_z_sols(:)))])
title('\Im(\partial_z\phi)')
colorbar

nexttile(2)
zztarg = reshape(phi_sols(1,:),size(xxgrid));
h2 = pcolor(xxgrid, yygrid, imag(zztarg),'EdgeColor','interp');
shading interp;
hold on
contour(xxgrid,yygrid,H,[4 4],'k','ShowText','off',"LabelFormat","%0.1f m")
clim(0.8*[min(real(phi_sols(:))), max(real(phi_sols(:)))])
title('\Im(\phi)')
colorbar

for kk = 1:1:ntimes
    nexttile(1)
    zztarg = reshape(phi_z_sols(kk,:),size(xxgrid));
    set(h1, 'CData', imag(zztarg));
    shading interp;

    nexttile(2)
    zztarg = reshape(phi_sols(kk,:),size(xxgrid));
    set(h2, 'CData', imag(zztarg)); 
    shading interp;

    title(t,sprintf('t = %.3f', ts(kk)));

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);