function [icoefs,dinds,pcoefs,H] = doline(X,Y,amp,width1,width2,w,eps)

    E = 7*10^9;
    nu = 0.33;
    H0 = 240;
    rhow = 1028;
    rhoi = 900;
    g = 9.8;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;

    R = sqrt(X.^2 + Y.^2);
    Hbar = amp.*(R < width1) + (amp + amp.*(width1 - R)/(width2 - width1)) .* (R < width2) .* (R > width1);
    H = H0 + Hbar;

    dHdR = (-amp/(width2 - width1)) .* (R < width2) .* (R > width1);
    dH2dR2 = 0;

    dRdX = X.*R.^(-1);
    dRdY = Y.*R.^(-1);
    d2RdX2 = R.^(-1) - X.^2.*R.^(-3);
    d2RdY2 = R.^(-1) - Y.^2.*R.^(-3);
    d2RdXdY = -X.*Y.*R.^(-3);

    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;

    Hx = dHdR.*dRdX;
    Hy = dHdR.*dRdY;

    Hxx = dH2dR2.*dRdX.^2 + dHdR.*d2RdX2;
    Hxy = dH2dR2.*dRdX.*dRdY + dHdR.*d2RdXdY;
    Hyy = dH2dR2.*dRdY.^2 + dHdR.*d2RdY2;

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    pcoefs = {a0/E,abar/E,b0/E,bbar/E,g0/E,gbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

    dinds = find(abs(Hbar / H0) > eps);

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