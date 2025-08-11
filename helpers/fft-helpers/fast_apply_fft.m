function v = fast_apply_fft(mu,kerns,coefs,iinds,jinds,N)

    if size(coefs,3) ~= size(kerns,3)
        error('coefs and kerns must have the same number of pages')
    end


    dinds_aug = sub2ind([N N], iinds, jinds);

    mu_aug = zeros(N);
    mu_aug(dinds_aug) = mu;
    mu_hat = fft2(mu_aug);

    G_mu = ifft2(kerns .* mu_hat);
    G_mu = reshape(G_mu,[N^2,1,size(kerns,3)]);
    G_mu = G_mu(dinds_aug,1,:);
    v = mu + sum(coefs .* G_mu,3);

end