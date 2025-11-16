function v = fast_apply_fft_nonlinear_flex(delta_mu_k,u_k,zk,a0,a1,kerns,coefs,iinds,jinds,N)

    if size(coefs,3) ~= size(kerns,3)
        error('coefs and kerns must have the same number of pages')
    end


    dinds_aug = sub2ind([N N], iinds, jinds);

    mu_aug = zeros(N);
    mu_aug(dinds_aug) = delta_mu_k;
    mu_hat = fft2(mu_aug);

    G_delta_mu = ifft2(kerns .* mu_hat);
    G_delta_mu = reshape(G_delta_mu,[N^2,1,size(kerns,3)]);
    G_delta_mu = G_delta_mu(dinds_aug,1,:);
    v = delta_mu_k - zk.^4.*coefs.*(a0 + 2*a1.*u_k.*conj(u_k)).*G_delta_mu ...
        -  zk.^4.*coefs.*a1.*u_k.*u_k.*conj(G_delta_mu);

end