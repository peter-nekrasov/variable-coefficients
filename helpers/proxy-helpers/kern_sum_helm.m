function out = kern_sum_helm(zk,src,targ)

V = targ.V;

kerns = helm2d.green_cell_helm(zk,src.r,targ.r);

G = kerns{1};

out = zk^2*V.*G ;

end