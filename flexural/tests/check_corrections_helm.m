%% Checking corrections for G

zk = 0.9;
targ = [2; 2];

dens = @(x,y) 10*x.*exp(-(x.^4+y.^4)/(10));
truev = integral2(@(x,y) dens(x,y).*helmgreenvalonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];
errs0 = hs*0;
errs1 = hs*0; 

gfunc = @(s,t) helmgreen1(zk,s,t);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_helm(h,zk);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

figure(1);clf

loglog(hs,errs0,'o-')
hold on

loglog(hs,errs1,'o-')
% hold on

loglog(hs,0.02*hs.^2,'--')
hold on

loglog(hs,0.0005*hs.^4,'--')
hold on

loglog(hs,0.0005*hs.^6,'--')
hold on

loglog(hs,0.0005*hs.^8,'--')
hold on

loglog(hs,errs1,'o-','DisplayName','$K_1, K_2$')
xlim([hs(end) hs(1)])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')


function val = helmgreenvalonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = helmgreen1(zk,src,targ);
    val = reshape(gf(:,:,1),size(x));
    

end
