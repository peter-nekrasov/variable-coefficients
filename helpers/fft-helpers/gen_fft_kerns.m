function fft_kerns = gen_fft_kerns(kerns,sz,ind)

    [zi,zj] = ind2sub(sz,ind);
    % fft_kerns = zeros([sqrt(size(kerns,1)) sqrt(size(kerns,1)) size(kerns,3)]);

    for ii = 1:size(kerns,3)
        GS = reshape(kerns(:,:,ii),sz);
        GS = circshift(GS,-[zi-1,zj-1]);
        fft_kerns(:,:,ii) = fft2(GS);
    end
    
end