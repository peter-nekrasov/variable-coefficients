function [field_sensors, field_domain] = forward_problem_flex_sd(X,Y,kerns,zk,dir,sensors,q,dx)

lmax = max(X(1,:));
N1 = size(X,1);
N2 = 2^ceil(log(ceil(4*lmax/dx)) / log(2));

dinds = find(abs(q(:)) > 1e-12 );
[iinds,jinds] = ind2sub(size(X),dinds);

theta = atan2(dir(2),dir(1));

coefs = -zk^4*q(dinds);
uinc = planewave(zk,[X(:) Y(:)].',theta);
rhs_vec = get_rhs(coefs,uinc,dinds);

sol = gmres(@(mu)  fast_apply_fft(mu,kerns,coefs,iinds,jinds,N2),rhs_vec,[],1e-10,size(rhs_vec,1));
% mu = zeros(size(xxgrid));
% mu(dinds) = sol;

usca = sol_eval_fft(sol,kerns,iinds,jinds,N1,N2);

utot = usca + uinc;
field_domain = reshape(utot,size(X));

field_sensors = dx^2 * (flexgreen(zk,[X(dinds) Y(dinds)].',sensors)*sol);

end