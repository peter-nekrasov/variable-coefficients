L = 1000; % length of grid
N1 = 400; % number of grid points

xs = linspace(-L/2,L/2,N1);
[xxgrid,yygrid] = meshgrid(xs);
h = xs(2) - xs(1);

[coefs,dinds,pcoefs,H,H0] = bump2(xxgrid,yygrid,5,50,8,1e-8);
[coefs2,dinds2,pcoefs2,H] = make_coefs(xxgrid,yygrid,H0,H-H0);

vecnorm(coefs(:) - coefs2(:)) ./ vecnorm(coefs(:))

vecnorm(dinds(:) - dinds2(:)) ./ vecnorm(dinds(:))

for ii = 1:12
    vecnorm(pcoefs{ii}(:) - pcoefs2{ii}(:)) 
end