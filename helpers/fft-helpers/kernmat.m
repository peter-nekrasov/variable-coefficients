function kerns = kernmat(src,targ,func,h,inds,corrs)
% Evaluates kernels and adds corrections to the source part based on inds

    kerns = func(src,targ);

    if nargin < 5
        correct = false;
    else
        correct = true;
    end

    if correct
        if (size(kerns,3) ~= numel(inds)) | (numel(inds) ~= numel(corrs))
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

    kerns = kerns*h*h;
    
    if correct
        for ii = 1:length(corrs)
            ind = inds{ii};
            corr = corrs{ii};
            kern = kerns(:,:,ii);
            for jj = 1:numel(corr(:,:,ii))
                kern((round(dx/h) == ind(jj,1)) & (round(dy/h) == ind(jj,2))) = ...
                    kern((round(dx/h) == ind(jj,1)) & (round(dy/h) == ind(jj,2))) + h*h*corr(jj,1);
            end
            kerns(:,:,ii) = kern;
        end
    end

end