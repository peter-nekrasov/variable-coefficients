figure(1); clf;
load("average_field_1.mat")

t = tiledlayout(1,2,"TileSpacing","compact");
nexttile
pcolor(xxgrid,yygrid,phizavg,'EdgeColor','none')
colorbar
title("|\phi_z|^2")

nexttile
pcolor(xxgrid,yygrid,lapphizavg,'EdgeColor','none')
colorbar
title("|\Delta\phi_z|^2")

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
