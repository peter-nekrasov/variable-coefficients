function [icoefs,dinds,pcoefs,H] = bumps(X,Y,xmin,xmax,amp,zk,eps)

    E = 7*10^9;
    nu = 0.33;
    H0 = 1;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;
    D = E*H0^3/(12*(1-nu^2));

    w = sqrt((D*zk^5 + rhow*g*zk) / (rhow*(H0*zk+1)));

    a0 = D;
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;

    [Xg,Yg] = meshgrid(linspace(xmin,xmax,64));
    sigma = (xmax-xmin)/(64);

    Xg = Xg(:);
    Yg = Yg(:);

    H = H0;
    Hx = 0;
    Hy = 0;
    Hxx = 0;
    Hxy = 0;
    Hyy = 0;

    for ii = 1:numel(Xg)

    Xo = Xg(ii);
    Yo = Yg(ii);

    Hbar = 2*(rand-0.5)*amp*exp(-((X-Xo).^2 + (Y-Yo).^2)/(sigma^2));

    H = H + Hbar;

    tpHx = -2*(X-Xo).*Hbar/sigma^2;
    tpHy = -2*(Y-Yo).*Hbar/sigma^2;

    tpHxx = - 2*Hbar/sigma^2 - 2*(X-Xo).*tpHx/sigma^2;
    tpHxy = - 2*(X-Xo).*tpHy/sigma^2;
    tpHyy = - 2*Hbar/sigma^2 - 2*(Y-Yo).*tpHy/sigma^2;

    Hx = Hx + tpHx;
    Hy = Hy + tpHy;

    Hxx = Hxx + tpHxx;
    Hxy = Hxy + tpHxy;
    Hyy = Hyy + tpHyy;

    end

    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    pcoefs = {a0/E,abar/E,b0/E,bbar/E,g0/E,gbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

    dinds = find(abs((H-H0) / H0) > eps);

    icoefs = zeros([size(dinds) 8]);
    icoefs(:,:,1) = (abar(dinds)*b0 - a0*bbar(dinds)) / a0;
    icoefs(:,:,2) = -(1-nu)*ayy(dinds);
    icoefs(:,:,3) = 2*(1-nu)*axy(dinds);
    icoefs(:,:,4) = -(1-nu)*axx(dinds);
    icoefs(:,:,5) = (axx(dinds)+ayy(dinds));
    icoefs(:,:,6) = 2*ax(dinds);
    icoefs(:,:,7) = 2*ay(dinds);
    icoefs(:,:,8) = -g0*abar(dinds)/a0;
    icoefs = icoefs*a0 ./ alpha(dinds);
    icoefs = icoefs/E;

end