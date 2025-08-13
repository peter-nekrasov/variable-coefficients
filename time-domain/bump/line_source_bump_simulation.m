% Time domain simulation of a line of sources

% Spatial parameters

bump_amp = -2.5;
bump_wid = 75;

srcloc = [-300;-300];
theta = pi/3;

L = 1000;
N1 = 200;

% Frequency/time parameters

t0 = 0;
w0 = 10;
sigma = 2;

gtilde = @(w) exp(1i*t0*(w-w0)).*exp(-(w-w0).^2/(2*sigma.^2));

nw = 80;
ws = linspace(w0 - 4*sigma,w0 + 4*sigma,nw);

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);
[src,targ,ind,sz,N2] = get_fft_grid(N1,L,1);

[~,dinds,pcoefs,H] = bump2(xxgrid,yygrid,bump_amp,bump_wid,ws(1),1e-8);
[iinds,jinds] = ind2sub(size(xxgrid),dinds);
a0 = pcoefs{1}; 
[inds,corrs] = get_correct_fg(h,a0);

fprintf('Number of points: %d \n',size(dinds,1))

ints_phi = zeros(nw, length(xxgrid(:)));
ints_phi_z = zeros(nw, length(xxgrid(:)));

for ii = 1:nw
    disp(ii)
    [coefs,~,pcoefs] = bump2(xxgrid,yygrid,bump_amp,bump_wid,ws(ii),1e-8); 

    a0 = pcoefs{1}; 
    b0 = pcoefs{3}; 
    g0 = pcoefs{5}; 

    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    zk = rts((imag(rts) == 0) & (real(rts) > 0));
    ejs = ejs/a0;

    % uincs = green1d(rts,ejs,srcloc,[xxgrid(:) yygrid(:)].',theta);
    uincs = fggreen(srcloc,[xxgrid(:) yygrid(:)].',rts,ejs);
    % cnst = max(kerns{1}(:));
    % kerns = cellfun(@(x) x/cnst,kerns,'UniformOutput',false);
    phizinc = uincs(:,:,1)/2;
    phiinc = uincs(:,:,end);
    rhs_vec = get_rhs(coefs,uincs,dinds);
    coefs(:,:,1:end-1) = 1/2*coefs(:,:,1:end-1);

    gfunc = @(s,t) fggreen(s,t,rts,ejs);
    kerns = kernmat(src,targ,gfunc,h,inds,corrs);
    kerns = gen_fft_kerns(kerns,sz,ind);
    
    evalkerns = kerns(:,:,[1 8]);

    % Solve with GMRES
    start = tic;
    sol = gmres(@(mu) fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-11,200);
    mu = zeros(size(xxgrid));
    mu(dinds) = sol;
    t1 = toc(start);
    fprintf('%5.2e s : time to solve\n',t1)
    
    usca = sol_eval_fft(sol,evalkerns,iinds,jinds,N1,N2);
    phizsca = usca(:,:,1)/2;
    phisca = usca(:,:,2);
        
    phitot = phisca + phiinc;
    phiztot = phizsca + phizinc;
    

    % figure(1); clf
    % tiledlayout(1,2,'TileSpacing','compact');
    % phitotplot = reshape(phitot,size(xxgrid));
    % phiztotplot = reshape(phiztot,size(xxgrid));
    % nexttile
    % pcolor(xxgrid,yygrid,real(phitotplot),'EdgeColor','interp')
    % colorbar
    % nexttile
    % pcolor(xxgrid,yygrid,real(phiztotplot),'EdgeColor','interp')
    % colorbar
    % title(['\omega = ',ws(ii)]);
    % drawnow
    utots = cat(3,reshape(phiztot,size(xxgrid)),reshape(phitot,size(xxgrid)));
    [abs_err, rel_err] = get_fin_diff_err(xxgrid,yygrid,utots,h,pcoefs,50,50,zk,dinds,'fg')

    ints_phi(ii,:) = phitot.*gtilde(ws(ii));
    ints_phi_z(ii,:) = phiztot.*gtilde(ws(ii));
end

% figure(1); clf
% zztarg = reshape(ints(1,:),size(xxgrid));
% pcolor(xxgrid,yygrid,real(zztarg),'EdgeColor','interp');

%%

ntimes = 500;
ts = linspace(-1,6,ntimes);
expmat = exp(-1i*ts(:)*ws);

phi_sols = expmat*ints_phi;
phi_z_sols = expmat*ints_phi_z;

%%

v = VideoWriter('line_source_bump_scattering.mp4', 'MPEG-4');  
v.FrameRate = 30;  % adjust as desired
open(v);

f = figure('Position',[1 1 1288 544]); clf
t = tiledlayout(1,2,'TileSpacing','compact');

nexttile(1)
zztarg = reshape(phi_z_sols(1,:),size(xxgrid));
h1 = pcolor(xxgrid, yygrid, real(zztarg),'EdgeColor','interp');
shading interp;
hold on
contour(xxgrid,yygrid,H,[3.5 4.5],'k','ShowText','off',"LabelFormat","%0.1f m")
clim(0.8*[min(real(phi_z_sols(:))), max(real(phi_z_sols(:)))])
title('\partial_z \phi')
colorbar

nexttile(2)
zztarg = reshape(phi_sols(1,:),size(xxgrid));
h2 = pcolor(xxgrid, yygrid, real(zztarg),'EdgeColor','interp');
shading interp;
hold on
contour(xxgrid,yygrid,H,[3.5 4.5],'k','ShowText','off',"LabelFormat","%0.1f m")
clim(0.8*[min(real(phi_sols(:))), max(real(phi_sols(:)))])
title('\phi')
colorbar

for kk = 1:1:ntimes
    nexttile(1)
    zztarg = reshape(phi_z_sols(kk,:),size(xxgrid));
    set(h1, 'CData', real(zztarg));
    shading interp;

    nexttile(2)
    zztarg = reshape(phi_sols(kk,:),size(xxgrid));
    set(h2, 'CData', real(zztarg)); 
    shading interp;

    title(t,sprintf('t = %.3f', ts(kk)));

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);