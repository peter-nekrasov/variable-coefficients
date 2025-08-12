function [err1,err2] = get_fin_diff_err(X,Y,utots,h,coefs,xloc,yloc,zk,dinds,eqtype)
% finite difference test for helmholtz equation
% checks error of utot at the closest point to (xloc,yloc)

    [~,ind] = min((X(:) - xloc(:)).^2 + (Y(:) - yloc(:)).^2);
    [ii, jj] = ind2sub(size(X),ind);

    if strcmpi(eqtype,'helm')

    % d2 - partial_{xx} (4th order)
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
    zk2V(dinds) = coefs(:,:,1);

    % Residual error of total solution 
    usub = utots(ii-4:ii+4,jj-4:jj+4);
    term1 = sum(lap.*usub,'all');
    term2 = (zk^2 + zk2V(ind))*utots(ii,jj);
    err1 = abs(term1+term2) ;
    err2 = err1 / max(abs([term1,term2])) ;

    elseif strcmpi(eqtype,'fg')

    % finite difference stencils
    d4 = [1	-4	6	-4	1].';
    d4 = d4/h^4; 
    
    d2 = [0 1	-2	1 0].';
    d2 = d2 / h^2;

    % d3 - partial_{xxx} (2nd order)
	d3 = [	-1/2	1	0	-1	1/2].' / h^3;	

    % d1 - partial_{x} (2nd order)
    d1 = [ 0	-1/2	0	1/2  0].' / h;

    lap = zeros(5);
    lap(3,:) = d2.';
    lap(:,3) = lap(:,3) + d2;
    
    bilap = zeros(5);
    bilap(:,3) = d4;
    bilap(3,:) = bilap(3,:) + d4.';
    bilap = bilap + 2*(d2*d2.');

    gradlapx = zeros(5);
    gradlapx(3,:) = d3.';
    gradlapx = gradlapx + d2*d1.';
    
    gradlapy = gradlapx.';

    hessxx = zeros(5);
    hessxx(3,:) = d2.' ;

    hessyy = zeros(5);
    hessyy(:,3) = d2 ;

    hessxy = d1*d1.' ;

    % Error of scattered part
    % phi_z_sub = phi_z(ii-4:ii+4,jj-4:jj+4);
    % first = alpha(ii,jj)*sum(bilap.*phi_z_sub,'all') ;
    % second = -beta(ii,jj)*phi_z(ii,jj);
    % third = gamma*phi(ii,jj);
    % rhs = k*bbar(ii,jj)*phiinc(ii,jj);
    % err = abs(first + second + third - rhs) 

    alpha = coefs{1} + coefs{2};
    beta = coefs{3} + coefs{4};
    g0 = coefs{5};
    alphax = coefs{7};
    alphay = coefs{8};
    alphaxx = coefs{9};
    alphaxy = coefs{10};
    alphayy = coefs{11};
    nu = coefs{end};
    alphalap = alphaxx + alphayy;

    phi_z = utots(:,:,1);
    phi = utots(:,:,2);
    
    % Residual error of total solution 
    phi_z_sub = phi_z(ii-2:ii+2,jj-2:jj+2);
    first = alpha(ii,jj).*sum(bilap.*phi_z_sub,'all') ;
    second = 2*alphax(ii,jj).*sum(gradlapx.*phi_z_sub,'all');
    third = 2*alphay(ii,jj).*sum(gradlapy.*phi_z_sub,'all');
    fourth = alphalap(ii,jj).*sum(lap.*phi_z_sub,'all');
    fifth = -(1-nu)*alphayy(ii,jj).*sum(hessxx.*phi_z_sub,'all');
    sixth = -(1-nu)*alphaxx(ii,jj).*sum(hessyy.*phi_z_sub,'all');
    seventh = 2*(1-nu)*alphaxy(ii,jj).*sum(hessxy.*phi_z_sub,'all');
    bterm = -beta(ii,jj).*phi_z(ii,jj);
    gterm = g0.*phi(ii,jj);
    err1 = abs(first + second + third + fourth + fifth + sixth + seventh + bterm + gterm);
    err2 = err1 / max(abs([first second third fourth fifth sixth seventh bterm gterm]));

    end
    
end