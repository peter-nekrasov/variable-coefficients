function v = fast_apply_fft(mu,kerns,coefs,corr,h,dinds,iinds,jinds,N)

    if length(coefs) ~= length(kerns)
        error('coefs and kerns must be cell arrays of equal length')
    end

    v = mu;

    for ii = 1:length(kerns)
    
    cs = coefs{ii};

    G_hat = kerns{ii};
    G_corr = corr{ii}(dinds,dinds);

    dinds_aug = sub2ind([N N], iinds, jinds);

    mu_aug = zeros(N);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    G_mu_aug = ifft2(G_hat.*mu_aug_hat);
    G_mu = G_mu_aug(dinds_aug);
    G_mu = G_mu(:) + G_corr*mu*h*h;

    v = v + cs.*G_mu;

    end

end