function pxy = get_proxy_ann(rmin,nrad,ntheta,rmax)
    if nargin < 2
        nrad = 6;
    end
    if nargin < 3
        ntheta = 31;
        ntheta = 51;
    end
    if nargin < 4
        rmax = 1.5;
        % rmax = 5;
    end

    thetas = linspace(0,2*pi,ntheta+1); thetas = thetas(1:end-1);
    % rads = linspace(rmin, rmin*rmax,nrad).';
    rads = 1./linspace(1/(rmax), 1/rmin,nrad).';
    
    pxy_x = rads .* cos(thetas);
    pxy_y = rads .* sin(thetas);
    pxy = [pxy_x(:).'; pxy_y(:).'];
end
