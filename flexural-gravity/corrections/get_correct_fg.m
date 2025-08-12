function [inds, corrs] = get_correct_fg(h,a0)
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

    A5 = [1 1 1 1 1;
        0 1 -1 0 0;
        0 0 0 1 -1;
        0 1 1 0 0;
        0 0 0 1 1];

    % kernmat(zi,zj) = kernmat(zi,zj) + tau(1);
    % kernmat(zi,zj+1) = kernmat(zi,zj+1) + tau(2);
    % kernmat(zi,zj-1) = kernmat(zi,zj-1) + tau(3);
    % kernmat(zi+1,zj) = kernmat(zi+1,zj) + tau(4);
    % kernmat(zi-1,zj) = kernmat(zi-1,zj) + tau(5);
    % kernmat(zi,zj+2) = kernmat(zi,zj+2) + tau(6);
    % kernmat(zi,zj-2) = kernmat(zi,zj-2) + tau(7);
    % kernmat(zi+2,zj) = kernmat(zi+2,zj) + tau(8);
    % kernmat(zi-2,zj) = kernmat(zi-2,zj) + tau(9);
    % kernmat(zi+1,zj+1) = kernmat(zi+1,zj+1) + tau(10);
    % kernmat(zi+1,zj-1) = kernmat(zi+1,zj-1) + tau(11);
    % kernmat(zi-1,zj+1) = kernmat(zi-1,zj+1) + tau(12);    
    % kernmat(zi-1,zj-1) = kernmat(zi-1,zj-1) + tau(13);

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

    i1 = [0 0; 1 0; -1 0; 0 1; 0 -1];
    i2 = [0 0; 1 0; -1 0; 0 1; 0 -1; 2 0; -2 0; 0 2; 0 -2; 1 1; -1 1; 1 -1; -1 -1];
    
    % 1/2 r^2 log r^2 correction
    
    % [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    % z1 = imag(z1)*1e12;
    % [~,~,z2] = epstein_zeta(-4+1i*10^-12,1,0,1,1,0,0) ;
    % z2 = imag(z2)*1e12;
    % [~,~,z3] = epstein_zeta(-4+1i*10^-12,1,0,1,0,1,0) ;
    % z3 = imag(z3)*1e12;
    % 
    % b = [2*z1; 0; 0; z2/2+z3/8; z2/2+z3/8];
    % b = b*h^4;
    % tau = A5 \ b;

    tau = [-0.127635425068666;0.007612114264598;0.007612114264598;0.007612114264598;0.007612114264598];
    tau = tau*h^4;
    c0 = 1/(4*pi*a0); % change this back to 1/(4 pi)

    valcor = c0*tau;

    % log(|r|^2) + 2 x^2 / r^2 + 1 
    
    % [z0] = epstein_zeta(0+1i*10^-12,1,0,1) ;
    % z0 = imag(z0)*1e12 ;
    % [~,z1] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    % z1 = imag(z1)*1e12;

    z0 =   -1.310532925911509;
    z1 =   -0.048593484005136;

    b = [z0 + log(h); 0; 0; z1; z1];
    tau0 = A5 \ b;

    %tau0 = [-1.906493138461182;-0.024296742002568;-0.024296742002568;-0.024296742002568;-0.024296742002568];
    tau0 = tau0*h^2;

    % [~,z1] = epstein_zeta(0+1i*10^-12,1,0,1,1,0,0) ;
    % z1 = imag(z1)*1e12;
    % [~,~,z2] = epstein_zeta(-2+1i*10^-12,1,0,1,0,1,0) ;
    % z2 = imag(z2)*1e12;
    % [~,~,z3] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    % z3 = imag(z3)*1e12;
    % 
    % b = [2*z1; 0; 0; 2*z3; 1/2*z2];
    %b = b*h^2;
    %tau1 = A5 \ b; 
    tau1 = [0.500000000000000;-0.038111783397667;-0.038111783397667;0.038111783397667;0.038111783397667];
    tau1 = tau1*h^2;

    hessxxcor = [c0*(2*tau0 + 2*tau1 ); zeros(8,1)];

    % 2*log(|r|) + 2 y^2 / r^2 + 1 

    %b = [2*z1; 0; 0; 1/2*z2; 2*z3];
    %b = b*h^2;
    %tau1 = A5 \ b; 
    tau1 = [0.500000000000000;0.038111783397667;0.038111783397667;-0.038111783397667;-0.038111783397667];
    tau1 = tau1*h^2;

    hessyycor = [c0*(2*tau0 + 2*tau1 ); zeros(8,1)];

    % 2*x*y / r^2

    % [~,~,z2] = epstein_zeta(-2+1i*10^-12,1,0,1,0,1,0) ;
    % z2 = imag(z2)*1e12;

    % z2 =    0.152447133590667;

    % b = [0; 0; 0; 0; 0; 1/2*z2; 0; 0; 0; 0; 0; 0; 0];
    %b = b*h^2;
    %tau = A13 \ b; 
    tau = [zeros(9,1); 0.0190558916988334; -0.0190558916988334; -0.0190558916988334; 0.0190558916988334];
    tau = tau*h^2;

    hessxycor = 2*c0*tau;

    % 4 x / r^2

    % [~,z1] = epstein_zeta(0+1i*10^-12,1,0,1,1,0,0) ;
    % z1 = imag(z1)*1e12;
    % [~,~,z2] = epstein_zeta(-2+1i*10^-12,1,0,1,0,1,0) ;
    % z2 = imag(z2)*1e12;
    % [~,~,z3] = epstein_zeta(-2+1i*10^-12,1,0,1,1,0,0) ;
    % z3 = imag(z3)*1e12;
    % 
    % b = [0; 2*z1; 0; 0; 0; 0; 0; 1/2*z2; 2*z3; 0; 0; 0; 0];
    %b = b*h;
    %tau = A13 \ b;
    tau = [0; 0.307925477734889; -0.307925477734889; 0; 0; -0.0480186305662778; ...
        0.0480186305662778; 0; 0; 0.0190558916988334; -0.0190558916988334; ...
        0.0190558916988334; -0.0190558916988334];
    tau = tau*h;

    gradlapxcor = 4*c0*tau;

    % 4 y / r^2

    %b = [0; 0; 2*z1; 0; 0; 0; 1/2*z2; 0; 0; 2*z3; 0; 0; 0];
    %b = b*h;
    % tau = A13 \ b;
    tau = [0; 0; 0; 0.307925477734889; -0.307925477734889; 0; 0; ...
        -0.0480186305662778; 0.0480186305662778; 0.0190558916988334; ...
        0.0190558916988334; -0.0190558916988334; -0.0190558916988334];
    tau = tau*h;

    gradlapycor = 4*c0*tau;

    % |r|^3 correction

    %[z0,~] = epstein_zeta_int(-3,1,0,1,0,0,0) ;
    %[~,z1] = epstein_zeta_int(-5,1,0,1,1,0,0) ;

    %b = [-z0; 0; 0; -2/5*z1; -2/5*z1];
    %b = b*h^5;
    %tau = A5 \ b;

    tau = [-0.045568816967685;0.004043633831185;0.004043633831185;0.004043633831185;0.004043633831185];
    tau = tau*h^5;
    c0 = 1/4/4/gamma(1+3/2)^2/a0; 

    phicor = 1/2*c0*tau;

    valcor = valcor / h^2;
    phicor = phicor / h^2;
    hessxxcor = hessxxcor / h^2;
    hessxycor = hessxycor / h^2;
    hessyycor = hessyycor / h^2;
    gradlapxcor = gradlapxcor / h^2;
    gradlapycor = gradlapycor / h^2;
    
    inds = {i1,i2,i2,i2,i2,i2,i2,i1};
    corrs = {valcor,hessxxcor,hessxycor,hessyycor,hessxxcor+hessyycor,gradlapxcor,gradlapycor,phicor};



end