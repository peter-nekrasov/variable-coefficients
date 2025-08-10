function rhs = get_rhs(coefs,uincs,dinds)


if (size(coefs) ~= size(uincs))
    error('coefs and uincs must be the same size')
end

rhs =  - sum(coefs.*uincs(dinds,:,:),3);

end