function [inds, corr] = get_correct_helm(h,zk)
% Getting corrections for the kernels in the Lippman-Schwinger eq
%
% input:
% - rts: float vector - roots of polynomial
% - ejs: float vector - coefficients in partial fraction expansion
% - h: float - grid spacing
%
% output: 
% - inds: cell array - indices of corrections
% - corrs: cell array - corresponding corrections
%
% format of cell arrays mirrors that of green.m

    inds = cell(1);
    corr = cell(1);

    i2 = [0 0; 1 0; -1 0; 0 1; 0 -1; 2 0; -2 0; 0 2; 0 -2; 1 1; -1 1; 1 -1; -1 -1];

    inds{1} = i2;

    A5 = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];

    A13 = [ones(1,13); ...
        0 1 -1 0 0 2 -2 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 2 -2 1 1 -1 -1; ...
        0 1 1 0 0 4 4 0 0 1 1 1 1; ...
        0 0 0 1 1 0 0 4 4 1 1 1 1; ...
        zeros(1,9) 1 -1 -1 1; ...
        zeros(1,9) 1 1 -1 -1; ...
        zeros(1,9) 1 -1 1 -1; ...
        0 1 -1 0 0 8 -8 0 0 1 -1 1 -1; ...
        0 0 0 1 -1 0 0 8 -8 1 1 -1 -1; ...
        0 1 1 0 0 16 16 0 0 1 1 1 1; ...
        zeros(1,9) ones(1,4); ...
        0 0 0 1 1 0 0 16 16 1 1 1 1];

    % correction to 1/(2*pi)*log(|r|) using 13 point stencil

    % [z0] = imag(epstein_zeta(0+1i*10^-12,1,0,1))*1e12 ;
    % [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    % z1 = imag(z1)*1e12;
    % [~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
    % z2 = imag(z2)*1e12;
    % [~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
    % z3 = imag(z3)*1e12;
    
    z0 =   -1.310532925911509;
    z1 =   -0.048593484005136;
    z2 =    0.057568617195860;
    z3 =   -0.108480640549869;

    b = [z0 + log(h); 0; 0; z1; z1; 0; 0; 0; 0; 0; 1/2*z2; 1/8*z3; 1/2*z2];

    tau0 = A13 \ b;
    tau0 = 1/(2*pi)*tau0;

    %  correction to -1/(8*pi)r^2 log r using 5 point stencil

    b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
    tau1 = A5 \ b;
    tau1 = tau1*h^2;
    tau1 = -zk^2/(8*pi)*tau1;
    
    tau0(1:5) = tau0(1:5) +  tau1;
    corr{1} = tau0 ;

end