function sol = sol_eval_fft(mu,kerns,iinds,jinds,N,num)
% Apply kernels for evaluation using FFT

    dinds_aug = sub2ind([num num], iinds, jinds);

    mu_aug = zeros(num);
    mu_aug(dinds_aug) = mu;
    mu_aug_hat = fft2(mu_aug);

    sol = zeros([N^2,1,size(kerns,3)]);

    for ii = 1:size(kerns,3)

    G_aug_hat = kerns(:,:,ii);

    G_mu_aug = ifft2(G_aug_hat.*mu_aug_hat);
    G_mu = G_mu_aug(1:N,1:N);
    sol(:,1,ii) = G_mu(:);

    end

end
