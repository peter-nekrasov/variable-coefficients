%% Checking corrections for G 

zk = 0.9;
targ = [2; 2];

dens = @(x,y) 10*x.*exp(-(x.^2+y.^2)/(10));
truev = integral2(@(x,y) dens(x,y).*greenvalonly(targ,x,y),-10,10,-10,10,"AbsTol",1e-18,"RelTol",1e-18);

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4];
errs0 = hs*0;
errs1 = hs*0;

gfunc = @(s,t) flexgreen2(zk,s,t);


for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,gfunc,h);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev) / abs(truev);

    [inds, corrs] = get_correct_flex(h);
    kern = kernmat(src,targ,gfunc,h,inds,corrs);
    val = kern(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev)

end

figure(1); clf

plot(log10(hs),log10(errs0),'o-')
hold on

plot(log10(hs),log10(errs1),'o-')
hold on

plot(log10(hs),log10(0.02*hs.^4),'--')
hold on

% loglog(hs,0.0005*hs.^4,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^8,'--')
% hold on

% loglog(hs,errs1,'o-','DisplayName','$K_1, K_2$')
xlim([log10(hs(end)) log10(hs(1))])

legend('Without corrections','With corrections','h^2','h^4','h^6','h^8','Location','southeast')


return 

%% Checking corrections for \partial_x (\Delta G) 
zk = 3;
targ = [2; 2];

dens = @(x,y) 10*x.*exp(-(x.^4+y.^4)/(10));
truev = integral2(@(x,y) dens(x,y).*gradlapxonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18)

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{3};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{3};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end


% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% % hold on
% 
% loglog(hs,0.02*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on

K1 = errs1;
loglog(hs,K1,'o-','DisplayName','$K_1, K_2$')
hold on





%% Checking corrections for G_{xx} 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) 10*exp(-(x.^4+y.^4)/(10));
truev = -0.361849809804212 + 1.456685781616957i;
%integral2(@(x,y) dens(x,y).*hessxxonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18)

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{2};
    % val = val(:,:,1);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{2};
    val = val(:,:,1);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end


% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 


K6 = errs1;
loglog(hs,K6,'o-','DisplayName','$K_3, K_4, K_5$')
hold on


%% Checking corrections for G_{xy} 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) 10*sin(x+y).*exp(-(x.^4+y.^4)/(10));
truev = 1.494565115435828 + 1.279414908871585i; %integral2(@(x,y) dens(x,y).*hessxyonly(targ,x,y,rts,ejs),-10,10,-10,10,'AbsTol',10E-18,'RelTol',10E-18)

hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4 0.025/8];
errs0 = hs*0;
errs1 = hs*0; 

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-5:h:5);
    src = [X(:).'; Y(:).'];
    % kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    % val = kern{2};
    % val = val(:,:,2);
    % d1 = dens(X,Y);
    % dint = sum(val(:).*d1(:),'all');
    % errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{2};
    val = val(:,:,2);
    d1 = dens(X,Y);
    dint = sum(val(:).*d1(:),'all');
    errs1(ii) = abs(dint - truev) / abs(truev);

end

% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^4,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
%
% legend('no correction', '13 pt correction', 'h^4', 'h^6','Location','northwest') 


K4 = errs1;
loglog(hs,K4,'o-','DisplayName','$K_6$')
hold on







%% Checking corrections for G_{yy} 

% addpath(genpath('..'))
% 
% gamma = -1;
% beta = 3;
% [rts,ejs] = find_roots(beta,gamma);
% targ = [2; 2];
% 
% dens = @(x,y) 100*exp(-(x.^2+y.^2)/(10));
% truev = -0.207870844101995 + 1.209128473565091i;
% 
% hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4];
% errs0 = hs*0;
% errs1 = hs*0; 
% 
% for ii = 1:numel(hs)
% 
%     h = hs(ii);
%     [X,Y] = meshgrid(-15:h:15);
%     src = [X(:).'; Y(:).'];
%     kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
%     val = kern{2};
%     val = val(:,:,3);
%     d1 = dens(X,Y);
%     dint = sum(val(:).*d1(:),'all');
%     errs0(ii) = abs(dint - truev);
% 
%     [inds, corrs] = get_correct(h,1);
%     kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
%     val = kern{2};
%     val = val(:,:,3);
%     d1 = dens(X,Y);
%     dint = sum(val(:).*d1(:),'all');
%     errs1(ii) = abs(dint - truev) / abs(truev);
% 
% end
% 
% 
% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.1*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% 
% K5 = errs1;
% 
% legend('no correction', '5 pt correction', 'h^2', 'h^6','Location','northwest') 



%% Checking corrections for \partial_y (\Delta G) 

% addpath(genpath('..'))
% 
% gamma = -1;
% beta = 3;
% [rts,ejs] = find_roots(beta,gamma);
% targ = [2; 2];
% 
% dens = @(x,y) y.*exp(-(x.^2+y.^2)/(10));
% truev =  -0.089213311677999 + 0.116604914243712i; % integral2(@(x,y) dens(x,y).*gradlapxonly(targ,x,y,rts,ejs),-15,15,-15,15,'AbsTol',0,'RelTol',10E-12)
% 
% hs = [2 1 0.5 0.2 0.1 0.05 0.025 0.025/2];
% errs0 = hs*0;
% errs1 = hs*0; 
% 
% for ii = 1:numel(hs)
% 
%     h = hs(ii);
%     [X,Y] = meshgrid(-15:h:15);
%     src = [X(:).'; Y(:).'];
%     kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
%     val = kern{3};
%     val = val(:,:,2);
%     d1 = dens(X,Y);
%     dint = sum(val(:).*d1(:),'all');
%     errs0(ii) = abs(dint - truev);
% 
%     [inds, corrs] = get_correct(h,1);
%     kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
%     val = kern{3};
%     val = val(:,:,2);
%     d1 = dens(X,Y);
%     dint = sum(val(:).*d1(:),'all');
%     errs1(ii) = abs(dint - truev) / abs(truev);
% 
% end
% 
% 
% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% % hold on
% 
% loglog(hs,0.02*hs.^2,'--')
% hold on
% 
% loglog(hs,0.0005*hs.^6,'--')
% hold on
% 
% K2 = errs1;
% 
% K3 = K1/2+K2/2;
% 
% legend('no correction', '13 pt correction', 'h^2', 'h^6','Location','northwest') 

%%  Checking corrections for phi 

addpath(genpath('..'))

gamma = -1;
beta = 3;
[rts,ejs] = find_roots(beta,gamma);
targ = [2; 2];

dens = @(x,y) 10*exp(-(x.^2+y.^2)/(5));
greenfac = @(x,y) phivalonly(targ,x,y,rts,ejs);
%truev = -0.062850332948632 + 0.069122770940408i; %  integral2(@(x,y) dens(x,y).*greenfac(x,y),-50,50,-50,50,"AbsTol",0,"RelTol",10E-14);
truev = -2.361687377424415 - 0.503053288734189i; % integral2(@(x,y) dens(x,y).*greenfac(x,y),-30,30,-30,30,"AbsTol",10E-18,"RelTol",10E-18)

hs = [2 1 0.5 0.2 0.1 0.05 0.025];
errs0 = hs*0;
errs1 = hs*0;

for ii = 1:numel(hs)

    h = hs(ii);
    [X,Y] = meshgrid(-30:h:30);
    src = [X(:).'; Y(:).'];
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h);
    val = kern{4};
    d1 = dens(X,Y);
    dint = sum(val.'.*d1(:));
    errs0(ii) = abs(dint - truev);

    [inds, corrs] = get_correct(h,1);
    kern = kernmat(src,targ,@(s,t) green(s,t,rts,ejs),h,inds,corrs);
    val = kern{4};
    dint = sum(val.'.*d1(:));
    errs1(ii) = abs(dint - truev) / abs(truev);

end

% loglog(hs,errs0,'o-')
% hold on
% 
% loglog(hs,errs1,'o-')
% hold on
% 
% loglog(hs,0.0005*hs.^5,'--')
% hold on
% 
% loglog(hs,0.0001*hs.^6,'--')
% hold on 
% 
% legend('no correction', '5 pt correction', 'h^5', 'h^6') 


K8 = errs1;
loglog(hs,K8,'o-','DisplayName','$K_8$')
hold on

%%

hs = [0.5 0.2 0.1 0.05 0.025 0.025/2 0.025/4];

loglog(hs,10*hs.^6,'k--','DisplayName','$h^6$')

legend('Interpreter','latex','Location','eastoutside')

ylabel('Relative error')
xlabel('h')

xlim([2E-3 3])
ylim([5E-16 2])

set(gca, 'FontSize',12)
fontname(gcf, 'CMU Serif')


%%

saveas(gcf,'intconv.fig','fig')
exportgraphics(gcf,'intconv.pdf','ContentType','vector')

function val = greenvalonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,1),size(x));
    
end


function val = gradlapxonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,6),size(x));
    
end

function val = gradlapyonly(targ,x,y,zk)

    src = [x(:).'; y(:).'];
    gf = flexgreen2(zk,src,targ);
    val = reshape(gf(:,:,7),size(x));
    
end