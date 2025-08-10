function sol = sol_eval_fft_sub(mu,kerns,corrs,h,dinds,iinds,jinds,N,num)
% Apply kernels for evaluation using FFT w/ corrections

    if size(kerns,3) ~= size(corrs,3)
        error('kerns and corrs must have the same number of pages')
    end

    dinds_aug = sub2ind([num num], iinds, jinds);

    mu_aug = zeros(num);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    mu0 = zeros(N^2,1);
    mu0(dinds) = mu;

    sol = zeros([N^2,1,size(kerns,3)]);

    for ii = size(kerns,3)

    G_aug_hat = kerns(:,:,ii);
    G_corr = corrs{ii};

    G_mu_aug = ifft2(G_aug_hat.*mu_aug_hat);
    G_mu = G_mu_aug(1:N,1:N);
    sol(:,:,ii) = G_mu(:) + G_corr*mu0*h*h;

    end

end
