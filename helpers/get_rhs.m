function rhs = get_rhs(coefs,uincs,dinds)

rhs = 0;

if (length(coefs) ~= length(uincs))
    error('coefs and uincs must be the same size')
end

for ii = 1:length(coefs)

cf = coefs{ii};
ui = uincs{ii};

rhs = rhs - cf.*ui(dinds);

end


end