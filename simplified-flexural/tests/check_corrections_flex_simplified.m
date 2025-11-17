%% Checking corrections for G

zk = 0.9;
targ = [2; 2];

dens = @(x,y) 10*x.*exp(-(x.^4+y.^4)/(10));
truev = integral2(@(x,y) dens(x,y).*flexgreenvalonly(targ,x,y,zk),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18);

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];
errs0 = hs*0;
errs1 = hs*0; 

gfunc = @(s,t) flexgreen(zk,s,t);

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct_flex_simplified(h);
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


function val = flexgreenvalonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen(zk,src,targ);
    val = reshape(gf(:,:,1),size(x));
    

end
