function [icoefs,dinds,pcoefs,H] = make_coefs(X,Y,H0,Hbar,type)

    w = 8;
    eps = 1e-8;
    E = 7*10^9;
    nu = 0.33;
    rhow = 1025;
    rhoi = 917;
    g = 9.8;

    a0 = E*H0^3/(12*(1-nu^2));
    b0 = (rhoi*H0*w^2 - rhow*g);
    g0 = -w^2*rhow;
    gbar = 0;
    H = H0 + Hbar;

    alpha = E*H.^3/(12*(1-nu^2));
    beta = (rhoi*H*w^2 - rhow*g);

    abar = alpha - a0;
    bbar = beta - b0;

    [Ny, Nx] = size(X);
    
    dx = X(1, 2) - X(1, 1);
    dy = Y(2, 1) - Y(1, 1);
    
    Lx = Nx * dx;
    Ly = Ny * dy;
    
    kx = (2*pi/Lx) * [0:ceil(Nx/2)-1, -floor(Nx/2):-1];
    ky = (2*pi/Ly) * [0:ceil(Ny/2)-1, -floor(Ny/2):-1];
    
    [KX, KY] = meshgrid(kx, ky);
    
    F_hat = fft2(Hbar);
    
    Hx = real(ifft2(1i * KX .* F_hat));
    Hy = real(ifft2(1i * KY .* F_hat));
    
    Hxx = real(ifft2(-KX.^2 .* F_hat));
    Hyy = real(ifft2(-KY.^2 .* F_hat));
    Hxy = real(ifft2(-KX .* KY .* F_hat));

    ax = 3*E*H.^2.*Hx/(12*(1-nu^2));
    ay = 3*E*H.^2.*Hy/(12*(1-nu^2));

    axx = 6*E*H.*Hx.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hxx/(12*(1-nu^2));
    axy = 6*E*H.*Hx.*Hy/(12*(1-nu^2)) + 3*E*H.^2.*Hxy/(12*(1-nu^2));
    ayy = 6*E*H.*Hy.^2/(12*(1-nu^2)) + 3*E*H.^2.*Hyy/(12*(1-nu^2));

    dinds = find(abs(Hbar / H0) > eps);

    if strcmpi(type,'fg') % flexural-gravity coefficients
    pcoefs = {a0/E,abar/E,b0/E,bbar/E,g0/E,gbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

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
    elseif strcmpi(type,'flex') % flexural coefficients

    pcoefs = {a0/E,abar/E,b0/E,bbar/E,ax/E,ay/E,axx/E,axy/E,ayy/E,nu}; 

    icoefs = zeros([size(dinds) 7]);
    % icoefs(:,:,1) = (abar(dinds)*b0 - a0*bbar(dinds)) / a0;
    icoefs(:,:,1) = (alpha(dinds)*b0 - a0*beta(dinds)) / a0;
    icoefs(:,:,2) = -(1-nu)*ayy(dinds);
    icoefs(:,:,3) = 2*(1-nu)*axy(dinds);
    icoefs(:,:,4) = -(1-nu)*axx(dinds);
    icoefs(:,:,5) = (axx(dinds)+ayy(dinds));
    icoefs(:,:,6) = 2*ax(dinds);
    icoefs(:,:,7) = 2*ay(dinds);
    icoefs = icoefs*a0 ./ alpha(dinds);
    icoefs = icoefs/E;

    end

end