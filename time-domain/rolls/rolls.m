function [icoefs,dinds,pcoefs,H] = rolls(X,Y,xl,xr,yl,yr,amp,width,freq,eps)

    E = 7*10^9;
    nu = 0.33;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;
    w = freq; % frequency

    freb = 0.35;

    H0 = -rhow/(rhoi-rhow)*freb ;

    [psi,psix,psiy,psixx,psixy,psiyy]  = bumpfunc(X,Y,xl,xr,yl,yr);

    S = (amp/2+amp/2*sin(2*pi*X/width)).*psi + freb ;
    B = rhoi/(rhoi-rhow)*S;
    H = S - B;

    Sx = amp*pi/width*cos(2*pi*X/width).*psi + (amp/2+amp/2*sin(2*pi*X/width)).*psix;
    Sy = (amp/2+amp/2*sin(2*pi*X/width)).*psiy;

    Sxx = -2*amp*pi^2/width^2*sin(2*pi*X/width).*psi + ...
        amp*pi/width*cos(2*pi*X/width).*psix + ...
        amp*pi/width*cos(2*pi*X/width).*psix + ...
        (amp/2+amp/2*sin(2*pi*X/width)).*psixx;
    Sxy = amp*pi/width*cos(2*pi*X/width).*psiy + (amp/2+amp/2*sin(2*pi*X/width)).*psixy;
    Syy = (amp/2+amp/2*sin(2*pi*X/width)).*psiyy;

    Hx = rhow/(rhow - rhoi).*Sx;
    Hy = rhow/(rhow - rhoi).*Sy;
    Hxx = rhow/(rhow - rhoi).*Sxx;
    Hxy = rhow/(rhow - rhoi).*Sxy;
    Hyy = rhow/(rhow - rhoi).*Syy;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;

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

    Hbar = H - H0;
    dinds = find(abs(Hbar / H0) >= eps);

    icoefs = zeros([size(dinds) 8]);
    icoefs(:,:,1) = (abar(dinds)*b0 - a0*bbar(dinds)) / a0;
    icoefs(:,:,2) = -(1-nu)*ayy(dinds);
    icoefs(:,:,3) = 2*(1-nu)*axy(dinds);
    icoefs(:,:,4) = -(1-nu)*axx(dinds);
    icoefs(:,:,5) = (axx(dinds)+ayy(dinds));
    icoefs(:,:,6) = 2*ax(dinds);
    icoefs(:,:,7) = 2*ay(dinds);
    icoefs(:,:,8) = -g0*abar(dinds)/a0;
    % icoefs = icoefs*a0 ./ alpha(dinds);
    icoefs = icoefs/E;

end