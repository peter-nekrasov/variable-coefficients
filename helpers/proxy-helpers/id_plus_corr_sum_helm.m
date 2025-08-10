function out = id_plus_corr_sum_helm(zk,coefs,corr,dinds,h)

V = coefs{1}(dinds);

G_corr = corr{1}(dinds,dinds)*h^2;

out = speye(size(G_corr)) + zk^2*V.*G_corr;

end