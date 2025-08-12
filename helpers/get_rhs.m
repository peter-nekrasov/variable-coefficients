function rhs = get_rhs(coefs,uincs,dinds)


if all(size(coefs) ~= size(uincs))
    error('coefs and uincs must be the same size')
end

if nargin < 3
    dinds = 1:size(coefs,1);
end

rhs =  - sum(coefs.*uincs(dinds,:,:),3);

end