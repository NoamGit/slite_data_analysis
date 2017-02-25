function [kernel, gfit] =  findKernelFromIm(im, scale, imSize, st_STA, params)
%findKernelFromIm(im, scale, imSize, st_STA)
% given an image and other params, it fits the best RF gaussian profile,
% finds the 8 closes pixels in the scaled image (to image original size)
% and extracts the temporal kernel
    
    
    if nargin < 5
        params = struct();
    end
    
    %parse inputs
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('showflag',false);
    p.addOptional('gfit',[]);
    p.parse(params);
    params = p.Results;
    gfit = params.gfit;
    
    % find gaussian fit 2 ( + tilt) to scaled image
    opts.tilted = true;
    height = imSize(1); width = imSize(2);
    [xi,yi] = meshgrid(1:width,1:height);
    im_scale = imresize(im((26:end-25),(1:end)),1/scale,'bicubic');
    gfit_0 = autoGaussianSurf(xi,yi,im_scale,opts);
    if(isempty(params.gfit))
        gfit = gfit_0;
        if(params.showflag)
            figure(3)
            set(gcf,'position',[ 2000 400 1100 450]);
            map = [[linspace(0,0,24)';linspace(0,1,26)'] [linspace(0,0,40)';linspace(0,1,10)'] linspace(0,0.6,50)'];
            colormap(map);
            subplot(1,2,1);imagesc(im_scale);
            subplot(1,2,2);imagesc(xi(:),yi(:),gfit.G);
        end
    end
    
    % find covariance matrix according to Sigma = RSSR^-1 see 
    % http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    R = [cos(gfit.angle),-sin(gfit.angle); sin(gfit.angle),cos(gfit.angle)];
    S = [gfit.sigmax 0 ;0 gfit.sigmay];
    Sigma = R * S * S / R;

    %finds 8 closest neighbors
    y = ceil((1:height));
    x = ceil((1:width));
    [X,Y] = meshgrid(x,y);
    coord = [X(:),Y(:)];
    X_hat = bsxfun(@minus,coord,[gfit.x0,gfit.y0]);
    X_hat2 = [X_hat(:,1).^2 repmat(prod(X_hat,2),1,2) X_hat(:,2).^2];
    inv_Sigma = inv(Sigma);
    mahal_dist = X_hat2 * inv_Sigma(:); % mahanabolis distance in vectorized way for data
    [~,I_mahal] = sort(mahal_dist,'ascend');
    pixel_i = coord(I_mahal(1:8),:);

    % find kernel on resized F(x,y,t)
    kernel = zeros(size(st_STA,3),1);
    for k = 1:size(st_STA,3)
    %     win = imcrop(rfm_blurr(:,:,k),[wx-round(width/2),wy-round(width/2),width,height]);
        win = st_STA(:,:,k);
        win_resize = imresize(win((26:end-25),(1:end)),1/scale,'bicubic');
        if(isempty(params.gfit))
            win_weighted = win_resize .* mapminmax_nc(gfit.G,0,1);
            win_lin = win_weighted(sub2ind(size(win_resize),pixel_i(:,2),pixel_i(:,1)));
        else
            % extracting the gaussian weights form the original gfit
            X_hat = bsxfun(@minus,coord,[gfit_0.x0,gfit_0.y0]);
            X_hat2 = [X_hat(:,1).^2 repmat(prod(X_hat,2),1,2) X_hat(:,2).^2];
            mahal_dist = X_hat2 * inv_Sigma(:); % mahanabolis distance in vectorized way for data
            [~,I_mahal] = sort(mahal_dist,'ascend');
            pixel_i_0 = coord(I_mahal(1:8),:);
            weight_map = mapminmax_nc(gfit_0.G,0,1);
            weight_vec = weight_map(sub2ind(size(win_resize),pixel_i_0(:,2),pixel_i_0(:,1)));
            
            win_weighted = win_resize;
            win_lin = win_weighted(sub2ind(size(win_resize),pixel_i(:,2),pixel_i(:,1))) .* weight_vec;
        end
        kernel(k) = mean(win_lin(:));
    end
    
end

