function err1 = get_fin_diff_err(X,Y,utot,h,coefs,xloc,yloc,zk,dinds,eqtype)
% finite difference test for helmholtz equation
% checks error of utot at the closest point to (xloc,yloc)

    [~,ind] = min((X(:) - xloc(:)).^2 + (Y(:) - yloc(:)).^2);
    [ii, jj] = ind2sub(size(X),ind);

    if strcmpi(eqtype,'helm')

    % d2 - partial_{xx} (8th order)
    d2 = zeros(9, 1);
    d2(1) = -1/560;
    d2(2) = 8/315;
    d2(3) = -1/5;
    d2(4) = 8/5;
    d2(5) = -205/72;
    d2(6) = 8/5;
    d2(7) = -1/5;
    d2(8) = 8/315;
    d2(9) = -1/560;

    lap = zeros(9);
    lap(5,:) = d2.';
    lap(:,5) = lap(:,5) + d2;
    lap = lap / h^2;

    zk2V = zeros(size(X));
    zk2V(dinds) = coefs{1};

    % Residual error of total solution 
    usub = utot(ii-4:ii+4,jj-4:jj+4);
    term1 = sum(lap.*usub,'all');
    term2 = (zk^2 + zk2V(ind))*utot(ii,jj);
    err1 = abs(term1+term2) / max(abs([term1,term2])) ;

    end
    
end