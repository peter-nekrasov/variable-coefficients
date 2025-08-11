function v = fast_apply_fft_sub(mu,kerns,coefs,idspmats,iinds,jinds,N)

    if size(coefs,3) ~= size(kerns,3)
        error('coefs and kerns must have the same number of pages')
    end

    v = idspmats*mu;

    dinds_aug = sub2ind([N N], iinds, jinds);

    mu_aug = zeros(N);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    for ii = 1:size(kerns,3)
    
    cs = coefs(:,:,ii);

    G_hat = kerns(:,:,ii);

    G_mu_aug = ifft2(G_hat.*mu_aug_hat);
    G_mu = G_mu_aug(dinds_aug);
    G_mu = G_mu(:) ;

    v = v + cs.*G_mu;

    end

end