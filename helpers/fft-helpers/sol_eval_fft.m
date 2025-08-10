function sol = sol_eval_fft(mu,kerns,corrs,h,dinds,iinds,jinds,N,num)

    if length(kerns) ~= length(corrs)
        error('kerns and corrs must be cell arrays of equal length')
    end

    dinds_aug = sub2ind([num num], iinds, jinds);

    mu_aug = zeros(num);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    mu0 = zeros(N^2,1);
    mu0(dinds) = mu;

    sol = cell(length(kerns));

    for ii = length(kerns)

    G_aug_hat = kerns{ii};
    G_corr = corrs{ii};

    G_mu_aug = ifft2(G_aug_hat.*mu_aug_hat);
    G_mu = G_mu_aug(1:N,1:N);
    sol{ii} = G_mu(:) + G_corr*mu0*h*h;

    end

end
