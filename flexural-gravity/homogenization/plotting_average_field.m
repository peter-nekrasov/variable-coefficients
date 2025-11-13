figure(1); clf;
load("average_field_1.mat")

t = tiledlayout(1,2,"TileSpacing","compact");
nexttile
pcolor(xxgrid,yygrid,log10(phizavg),'EdgeColor','none')
colorbar
title("log_{10}|\phi_z|^2")

nexttile
pcolor(xxgrid,yygrid,log10(lapphizavg),'EdgeColor','none')
colorbar
title("log_{10}|\Delta\phi_z|^2")

%%

figure(2); clf;
load("random_field_3671.mat")

t = tiledlayout(1,3,"TileSpacing","compact");
nexttile
pcolor(xxgrid,yygrid,abs(H),'EdgeColor','none')
colorbar
title("H")

nexttile
pcolor(xxgrid,yygrid,abs(phiztot),'EdgeColor','none')
colorbar
title("|\phi_z|")

nexttile
pcolor(xxgrid,yygrid,abs(lapphiztot),'EdgeColor','none')
colorbar
title("|\Delta\phi_z|")


%%

inds = (yygrid < 1) & (yygrid > 0) ; %& (xxgrid > 0);

xx = xxgrid(inds); phizslice = phizavg(inds);
plot(xx, (phizslice));
