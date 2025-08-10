function C = get_sparse_corr_mat(sz,inds,corrs)
% Takes corrections, indices and forms sparse matrices of size prod(sz) 
% by prod(sz). Source and target must both have size sz.
% [X, Y] = meshgrid(1:sz(1), 1:sz(2)); srcind = [X(:), Y(:)];


    % [k1s, k2s] = meshgrid(1:srcsz(1),1:srcsz(2));
    % [k1t, k2t] = meshgrid(1:targsz(1),1:targsz(2));
    % 
    % ks = sub2ind(srcsz, k1s(:), k2s(:));
    % kt = sub2ind(targsz, k1t(:), k2t(:));
    %
    % srcsz = size(srcx);
    % targsz = size(targx);
    % 
    % [k1s, k2s] = meshgrid(1:srcsz(1), 1:srcsz(2));
    % [k1t, k2t] = meshgrid(1:targsz(1), 1:targsz(2));
    % 
    % ks = sub2ind(size(k1s), k1s(:), k2s(:));
    % kt = sub2ind(size(k1t), k1t(:), k2t(:));

    C = cell(length(corrs),1);

    nx = sz(2);
    ny = sz(1);

    [k1t, k2t] = meshgrid(1:nx, 1:ny);

    for ii = 1:numel(corrs)
        corr = corrs{ii};
        ind = inds{ii};

        kt = repmat((1:(nx*ny)).',1,numel(corr));
        kt = kt(:);

        corr = repmat(corr.',nx*ny,1);
        corr = corr(:);

        k1s = k1t(:) - ind(:,1).'; % check 
        k2s = k2t(:) - ind(:,2).';

        k1s = k1s(:);
        k2s = k2s(:);

        keep = find((abs(k1s - (nx+1)/2) <= (nx-1)/2).*(abs(k2s - (ny+1)/2) <= (ny-1)/2 ));

        k1s = k1s(keep);
        k2s = k2s(keep);
        kt = kt(keep);
        corr = corr(keep);

        ks = sub2ind(sz, k2s, k1s); % k2 is row, k1 is column

        cspar = sparse(kt,ks,corr);

        C{ii} = cspar;
    end

end