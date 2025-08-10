function [coefs,dinds] = bump(X,Y,amp,width,ncs,eps)

    if nargin < 5
        ncs = 6;
    end

    V = amp*exp(-(X(:).^2 + Y(:).^2)/(2*width^2));

    Vx = - X(:).*V/width^2;
    Vy = - Y(:).*V/width^2;

    Vxx = - V/width^2 - X(:).*Vx/width^2;
    Vxy = - X(:).*Vy/width^2;
    Vyy = - V/width^2 - Y(:).*Vy/width^2;

    dinds = find(abs(V) > eps );

    coefs = cat(3,V(dinds),Vx(dinds),Vy(dinds),Vxx(dinds),Vxy(dinds),Vyy(dinds));
    coefs = coefs(:,:,1:ncs);

end