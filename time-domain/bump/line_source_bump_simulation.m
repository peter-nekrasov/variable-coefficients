% Time domain simulation of a line of sources

% Spatial parameters

bump_amp = -4;
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

[coefs,dinds,pcoefs] = bump2(xxgrid,yygrid,bump_amp,bump_wid,ws(1),1e-8);
[iinds,jinds] = ind2sub(size(xxgrid),dinds);
a0 = pcoefs{1}; 
[inds,corrs] = get_correct_fg(h,a0);

fprintf('Number of points: %d \n',size(dinds,1))

ints_phi = zeros(nw, length(xxgrid(:)));
ints_phi_z = zeros(nw, length(xxgrid(:)));

for ii = 1:nw
    [coefs,~,pcoefs] = bump2(xxgrid,yygrid,bump_amp,bump_wid,ws(ii),1e-8); 

    a0 = pcoefs{1}; 
    b0 = pcoefs{3}; 
    g0 = pcoefs{5}; 

    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    zk = rts((imag(rts) == 0) & (real(rts) > 0));
    ejs = ejs/a0;

    uincs = green1d(rts,ejs,srcloc,[xxgrid(dinds) yygrid(dinds)].',theta);
    % cnst = max(kerns{1}(:));
    % kerns = cellfun(@(x) x/cnst,kerns,'UniformOutput',false);
    phizinc = srcval(:,:,1);
    phiinc = srcval(:,:,end);
    rhs_vec = get_rhs(coefs,uincs);
    coefs = coefs / 2;

    gfunc = @(s,t) fggreen(s,t,rts,ejs);
    kerns = kernmat(src,targ,gfunc,h);
    kerns = gen_fft_kerns(kerns,sz,ind);
    
    evalkerns = kerns(:,:,[1 8]);

    % Solve with GMRES
    start = tic;
    mu = gmres(@(mu) fast_apply_fft(mu,kerns,coefs),rhs_vec,[],1e-5,200);
    mu = reshape(mu, size(xxtarg));
    t1 = toc(start);
    fprintf('%5.2e s : time to solve\n',t1)
    
    [phi, phi_z] = sol_eval_fft(mu,evalkerns);
    
    phi_tot = phi(:) + phiinc;
    phi_z_tot = phi_z(:) + phizinc;
    % % 
    % % figure(1); clf
    % % zztarg = reshape(mu,size(xxtarg));
    % % pcolor(xxtarg,yytarg,real(zztarg),'EdgeColor','interp')
    % % colorbar

    ints_phi(ii,:) = phi_tot.*gtilde(ws(ii));
    ints_phi_z(ii,:) = phi_z_tot.*gtilde(ws(ii));
end

% figure(1); clf
% zztarg = reshape(ints(1,:),size(xxtarg));
% pcolor(xxtarg,yytarg,real(zztarg),'EdgeColor','interp');


%%
ntimes = 300;
ts = linspace(-2,5,ntimes);
expmat = exp(-1i*ts(:)*ws);

phi_sols = expmat*ints_phi;
phi_z_sols = expmat*ints_phi_z;

%%

v = VideoWriter('line_source_bump_scattering.mp4', 'MPEG-4');  
v.FrameRate = 25;  % adjust as desired
open(v);

f = figure('Position',[1 1 1288 544]); clf
t = tiledlayout(1,2,'TileSpacing','compact');

nexttile(1)
zztarg = reshape(phi_z_sols(1,:),size(xxtarg));
h1 = pcolor(xxtarg, yytarg, real(zztarg),'EdgeColor','interp');
shading interp;
hold on
contour(xxtarg,yytarg,H,[2 3 4],'k','ShowText','off',"LabelFormat","%0.1f m")
clim(0.5*[min(real(phi_z_sols(:))), max(real(phi_z_sols(:)))])
title('\phi_z')
colorbar

nexttile(2)
zztarg = reshape(phi_sols(1,:),size(xxtarg));
h2 = pcolor(xxtarg, yytarg, real(zztarg),'EdgeColor','interp');
shading interp;
hold on
contour(xxtarg,yytarg,H,[2 3 4],'k','ShowText','off',"LabelFormat","%0.1f m")
clim(0.5*[min(real(phi_sols(:))), max(real(phi_sols(:)))])
title('\phi')
colorbar

for kk = 1:ntimes
    nexttile(1)
    zztarg = reshape(phi_z_sols(kk,:),size(xxtarg));
    set(h1, 'CData', real(zztarg));
    shading interp;

    nexttile(2)
    zztarg = reshape(phi_sols(kk,:),size(xxtarg));
    set(h2, 'CData', real(zztarg)); 
    shading interp;

    title(t,sprintf('t = %.3f', ts(kk)));

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);


function val = get_green(w,x0,y0,xxtarg,yytarg)

    [coefs, H] = bump2(X,Y,2,50,w); % remove gbar from coefs vector

    a0 = coefs{1}; 
    b0 = coefs{3}; 
    g0 = coefs{5}; 
    
    % Finding positive real roots
    [rts,ejs] = find_roots(b0 / a0, g0 / a0);
    k = rts((imag(rts) == 0) & (real(rts) > 0));
    ejs = ejs/a0;

    src = [x0;y0];
    targ = [xxtarg(:) yytarg(:)].';

    val = green(src,targ,rts,ejs);
    val = val{1};

end