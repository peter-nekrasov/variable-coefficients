function [icoefs,dinds,pcoefs,rhoi] = bump_variable_rho(X,Y,amp,width,w,eps)

    E = 7*10^9;
    nu = 0.33;
    H0 = 5;
    rhow = 1025;
    rhoi0 = 917;
    rhoibar = amp*exp(-(X.^2 + Y.^2)/(2*width^2));

    rhoi = rhoi0 + rhoibar;

    g = 9.8;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi0*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;
    % Hbar = amp*exp(-(X.^2 + Y.^2)/(2*width^2));
    Hbar = 0;
    H = H0 + Hbar;

    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;

    % Hx = -X.*Hbar/width^2;
    % Hy = -Y.*Hbar/width^2;
    % 
    % Hxx = - Hbar/width^2 - X.*Hx/width^2;
    % Hxy = - X.*Hy/width^2;
    % Hyy = - Hbar/width^2 - Y.*Hy/width^2;
    % 
    % ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    % ay = 3*E*H.^2.*Hy/(12*(1-nu^2));
    % 
    % axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    % axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    % ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    Hx = 0;
    Hy = 0;

    Hxx = 0;
    Hxy = 0;
    Hyy = 0;

    ax = 0;
    ay = 0;

    axx = 0;
    axy = 0;
    ayy = 0;

    pcoefs = {a0/E,abar/E,b0/E,bbar/E,g0/E,gbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

    dinds = find(abs(bbar / b0) > eps);

    icoefs = zeros([size(dinds) 8]);
    icoefs(:,:,1) = (abar*b0 - a0*bbar(dinds)) / a0;
    icoefs(:,:,2) = -(1-nu)*ayy;
    icoefs(:,:,3) = 2*(1-nu)*axy;
    icoefs(:,:,4) = -(1-nu)*axx;
    icoefs(:,:,5) = (axx+ayy);
    icoefs(:,:,6) = 2*ax;
    icoefs(:,:,7) = 2*ay;
    icoefs(:,:,8) = -g0*abar/a0;
    icoefs = icoefs*a0 ./ alpha;
    icoefs = icoefs/E;

end