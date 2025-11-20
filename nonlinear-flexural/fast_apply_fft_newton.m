function v = fast_apply_fft_newton(delta_mu_k,kerns,coefs,iinds,jinds,N)

    % if size(coefs,3) ~= size(kerns,3)
    %     error('coefs and kerns must have the same number of pages')
    % end

    dinds_aug = sub2ind([N N], iinds, jinds);

    delta_mu1 = delta_mu_k(1:numel(iinds));
    delta_mu2 = delta_mu_k(numel(iinds)+1:end);

    mu_aug1 = zeros(N);
    mu_aug1(dinds_aug) = delta_mu1;
    mu_hat1 = fft2(mu_aug1);

    mu_aug2 = zeros(N);
    mu_aug2(dinds_aug) = delta_mu2;
    mu_hat2 = fft2(mu_aug2);

    G_delta_mu1 = ifft2(kerns .* mu_hat1);
    G_delta_mu1 = reshape(G_delta_mu1,[N^2,1,size(kerns,3)]);
    G_delta_mu1 = G_delta_mu1(dinds_aug,1,:);

    G_delta_mu2 = ifft2(conj(fliplr(flipud(kerns))) .* mu_hat2);
    G_delta_mu2 = reshape(G_delta_mu2,[N^2,1,size(kerns,3)]);
    G_delta_mu2 = G_delta_mu2(dinds_aug,1,:);

    v1 = delta_mu1 + coefs(:,:,1).*G_delta_mu1 + coefs(:,:,2).*G_delta_mu2;
    v2 = delta_mu2 + coefs(:,:,1).*G_delta_mu2 + coefs(:,:,2).*G_delta_mu1;

    v = [v1; v2];

end