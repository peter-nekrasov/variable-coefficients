function [err1,err2] = get_fin_diff_err(X,Y,utots,h,coefs,xloc,yloc,zk,dinds,eqtype)
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
    zk2V(dinds) = coefs(:,:,1);

    % Residual error of total solution 
    usub = utots(ii-4:ii+4,jj-4:jj+4);
    term1 = sum(lap.*usub,'all');
    term2 = (zk^2 + zk2V(ind))*utots(ii,jj);
    err1 = abs(term1+term2) / max(abs([term1,term2])) ;

    elseif strcmpi(eqtype,'fg')

    % finite difference stencils 
    % d4 - partial_{xxxx} (6th order) % improve order?
    d4 = zeros(9,1);
    d4(1) = 7/240;	
    d4(2) = -2/5;
    d4(3) = 169/60;
    d4(4) = -122/15;
    d4(5) = 91/8;
    d4(6) = -122/15;
    d4(7) = 169/60;
    d4(8) = -2/5;
    d4(9) = 7/240;
    
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

    % d3 - partial_{xxx} (6th order)
	d3 = [-7/240; 3/10; -169/120; 61/30; 0; -61/30; 169/120; -3/10; 7/240];	

    % d1 - partial_{x} (8th order)
    d1 = [1/280; -4/105; 1/5; -4/5; 0; 4/5; -1/5; 4/105; -1/280];

    lap = zeros(9);
    lap(5,:) = d2.';
    lap(:,5) = lap(:,5) + d2;
    lap = lap / h^2;
    
    bilap = zeros(9);
    bilap(:,5) = d4;
    bilap(5,:) = bilap(5,:) + d4.';
    bilap = bilap + 2*(d2*d2.');
    bilap = bilap / h^4;

    gradlapx = zeros(9);
    gradlapx(5,:) = d3.';
    gradlapx = gradlapx + d2*d1.';
    gradlapx = gradlapx / h^3;
    
    gradlapy = gradlapx.';

    hessxx = zeros(9);
    hessxx(5,:) = d2.' / h^2;

    hessyy = zeros(9);
    hessyy(:,5) = d2 / h^2;

    hessxy = d1*d1.' / h^2;

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
    phi_z_sub = phi_z(ii-4:ii+4,jj-4:jj+4);
    first = alpha(ii,jj).*sum(bilap.*phi_z_sub,'all') ;
    second = 2*alphax(ii,jj).*sum(gradlapx.*phi_z_sub,'all');
    third = 2*alphay(ii,jj).*sum(gradlapy.*phi_z_sub,'all');
    fourth = alphalap(ii,jj).*sum(lap.*phi_z_sub,'all');
    fifth = -(1-nu)*alphayy(ii,jj).*sum(hessxx.*phi_z_sub,'all');
    sixth = -(1-nu)*alphaxx(ii,jj).*sum(hessyy.*phi_z_sub,'all');
    seventh = 2*(1-nu)*alphaxy(ii,jj).*sum(hessxy.*phi_z_sub,'all');
    bterm = -beta(ii,jj).*phi_z(ii,jj);
    gterm = g0.*phi(ii,jj);
    err1 = abs(first + second + third + fourth + fifth + sixth + seventh + bterm + gterm) ;
    err2 = err1 / max(abs([first second third fourth fifth sixth seventh bterm gterm]));

    elseif strcmpi(eqtype,'flex')

    % d4 - partial_{xxxx} (6th order) % improve order?
    d4 = zeros(9,1);
    d4(1) = 7/240;	
    d4(2) = -2/5;
    d4(3) = 169/60;
    d4(4) = -122/15;
    d4(5) = 91/8;
    d4(6) = -122/15;
    d4(7) = 169/60;
    d4(8) = -2/5;
    d4(9) = 7/240;
    
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
    
    bilap = zeros(9);
    bilap(:,5) = d4;
    bilap(5,:) = bilap(5,:) + d4.';
    bilap = bilap + 2*(d2*d2.');
    bilap = bilap / h^4;

    minzk4V = zeros(size(X));
    minzk4V(dinds) = coefs(:,:,1);

    % Residual error of total solution 
    usub = utots(ii-4:ii+4,jj-4:jj+4);
    term1 = sum(bilap.*usub,'all');
    term2 = (-zk^4 + minzk4V(ind))*utots(ii,jj);
    err1 = abs(term1+term2)  ;
    err2 = err1 / max(abs([term1 term2]));

    end

    
end