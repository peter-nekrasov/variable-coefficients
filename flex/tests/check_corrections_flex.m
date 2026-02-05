%% Checking corrections for G

zk = 0.9;
targ = [2; 2];

dens = @(x,y) 10*x.*exp(-(x.^4+y.^4)/(10));
truev = integral2(@(x,y) dens(x,y).*flexgreenvalonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

hs = [1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4];
errs0 = hs*0;
errs1 = hs*0; 

gfunc = @(s,t) flexgreen2(zk,s,t);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(1);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')


return 

%% Checking corrections for \partial_x (\Delta G) 

truev = integral2(@(x,y) dens(x,y).*gradlapxonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-14,'RelTol',10E-14);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,6);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,6);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(2);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')

%% Checking corrections for \partial_y (\Delta G) 

truev = integral2(@(x,y) dens(x,y).*gradlapyonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-14,'RelTol',10E-14);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,7);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,7);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(2);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')


%% Checking corrections for \partial_{xy} G 

truev = integral2(@(x,y) dens(x,y).*hessxyonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-14,'RelTol',10E-14);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,3);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,3);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(2);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')


%% Checking corrections for \partial_{xx} G 

truev = integral2(@(x,y) dens(x,y).*hessxxonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-14,'RelTol',10E-14);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(2);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')

%% Checking corrections for \partial_{yy} G 

truev = integral2(@(x,y) dens(x,y).*hessyyonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-14,'RelTol',10E-14);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,4);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,4);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(2);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')


%% Checking corrections for Laplacian of G 

truev = integral2(@(x,y) dens(x,y).*laponly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-14,'RelTol',10E-14);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,2) + kern(:,:,4);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,2) + kern(:,:,4);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(2);clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
% hold on

plot(log10(hs),log10(0.02*hs.^2),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^4),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^6),'--')
hold on

plot(log10(hs),log10(0.0005*hs.^8),'--')
hold on

plot(log10(hs),log10(errs1),'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')



function val = flexgreenvalonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,1),size(x));
    
end


function val = gradlapxonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,6),size(x));
    
end

function val = gradlapyonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,7),size(x));
    
end

function val = hessxxonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,2),size(x));
    
end

function val = hessxyonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,3),size(x));
    
end

function val = hessyyonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,4),size(x));
    
end

function val = laponly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,5),size(x));
    
end