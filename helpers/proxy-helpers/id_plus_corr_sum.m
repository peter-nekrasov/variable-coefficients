function out = id_plus_corr_sum(coefs,corr,dinds,h)

    if size(coefs,3) ~= length(corr)
        error('coefs and corr must be cell arrays of equal length')
    end

    out = speye(size(dinds,1));

    for ii = 1:length(corr)
    V = coefs(:,:,ii);
    G_corr = corr{ii}(dinds,dinds)*h^2;
    
    out = out + V.*G_corr;
    end

end