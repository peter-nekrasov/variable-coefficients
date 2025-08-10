function M = gen_fft_kerns(kerns,sz,ind)

    [zi,zj] = ind2sub(sz,ind);
    M = cell(1,numel(kerns));

    for ii = 1:numel(kerns)
        kern = kerns{ii};
        sz3 = size(kern,3);
        fftkern = zeros([sz, sz3]);
        for kk = 1:sz3
            GS = reshape(kern(:,:,kk),sz);
            GS = circshift(GS,-[zi-1,zj-1]);
            fftkern(:,:,kk) = fft2(GS);
        end
        M{ii} = fftkern;
    end
    
end