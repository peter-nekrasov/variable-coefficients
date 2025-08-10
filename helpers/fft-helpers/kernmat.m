function M = kernmat(src,targ,func,h,inds,corrs)
% Evaluates kernels and adds corrections to the source part based on inds

    kerns = func(src,targ);

    if nargin < 5
        correct = false;
    else
        correct = true;
    end

    if correct
        if (numel(kerns) ~= numel(inds)) | (numel(inds) ~= numel(corrs))
            error('output of func must be the same length as inds and corrs')
        end
    end

    [~,ns] = size(src);
    [~,nt] = size(targ);
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    xt = repmat(targ(1,:).',1,ns);
    yt = repmat(targ(2,:).',1,ns);
    
    dx = xt-xs;
    dy = yt-ys;
    
    dx2 = dx.*dx;
    dy2 = dy.*dy;
    
    r2 = dx2 + dy2;
    r = sqrt(r2);

    M = cell(1,numel(kerns));

    for ii = 1:numel(kerns)
        kern = kerns{ii};
        kern = kern*h*h;
        if correct
            ind = inds{ii};
            corr = corrs{ii};
            sz3 = size(corr,3);
            for kk = 1:sz3
                tmpker = kern(:,:,kk);
                tmpcor = corr(:,:,kk);
                for jj = 1:numel(corr(:,:,1))
                    tmpker((round(dx/h) == ind(jj,1)) & (round(dy/h) == ind(jj,2))) = ...
                        tmpker((round(dx/h) == ind(jj,1)) & (round(dy/h) == ind(jj,2))) + h*h*tmpcor(jj,1);
                end
                kern(:,:,kk) = tmpker;
            end
        end
        M{ii} = kern;
    end

end