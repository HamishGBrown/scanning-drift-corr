def makeImage(images,scanOr,scanDir indLines = None,indImages = None):
               

    """ This function generates a resampled scanning probe image with dimensions
    # of imageSize, from a an array of N scan lines given in scaneLines,
    # (lines specified as image rows), from an array of Nx2 origins in scanOr.
    # scanDir is a 2 element vector specifying the direction of the scan.
    # All arrays are stored inside struct sMerge.  ind specified update index.
    # indLines is a vector of binary values specifying which lines to include."""
    
    #Get number of images and image shape
    nimages,nopiy,nopix = images.shape
    
    # Boolean arrays indLine and indImages dictate which scan lines
    # and images are to be included in making the final image.
    # By default all images and lines are included
    if indLines  is None: 
        indLines  = np.ones((images.shape[2],),dtype=np.bool))
    if indImages is None: 
        indImages = np.ones((images.shape[0],),dtype=np.bool))
    
    from numpy import newaxis as na
    # Calculate pixel values with drift
    t = np.arange(images.shape[1])
    xInd = scanOr[indLines,indImages,1] + (t[na,:]*scanDir[indImage,1:])[na,:]
    yInd = scanOr[indLines,indImages,0] + (t[na,:]*scanDir[indImage,:1])[na,:]

    # Prevent pixels from leaving image boundaries
    xInd = np.clip(xInd,0,nopix-1)
    yInd = np.clip(yInd,0,nopiy-1)

    # Convert to bilinear interpolants and weights
    xIndF = np.floor(xInd)
    yIndF = np.floor(yInd)
    xAll = [xIndF xIndF+1 xIndF xIndF+1];
    yAll = [yIndF yIndF yIndF+1 yIndF+1];
    dx = xInd-xIndF;
    dy = yInd-yIndF;
    w = [(1-dx).*(1-dy) dx.*(1-dy) (1-dx).*dy dx.*dy];
    indAll = sub2ind(sMerge.imageSize,xAll,yAll);

    % Generate image
    sL = sMerge.scanLines(indLines,:,indImage);
    sig = reshape(accumarray(indAll(:),[ ...
        w(:,1).*sL(:);
        w(:,2).*sL(:);
        w(:,3).*sL(:);
        w(:,4).*sL(:)],...
        [prod(sMerge.imageSize) 1]),sMerge.imageSize);
    count = reshape(accumarray(indAll(:),[ ...
        w(:,1);w(:,2);w(:,3);w(:,4)],...
        [prod(sMerge.imageSize) 1]),sMerge.imageSize);

    % Apply KDE
    r = max(ceil(sMerge.KDEsigma*3),5);
    sm = fspecial('gaussian',2*r+1,sMerge.KDEsigma);
    sm = sm / sum(sm(:));
    sig = conv2(sig,sm,'same');
    count = conv2(count,sm,'same');
    sub = count > 0;
    sig(sub) = sig(sub) ./ count(sub);
    sMerge.imageTransform(:,:,indImage) = sig;

    % Estimate sampling density
    bound = count == 0;
    bound([1 end],:) = true;
    bound(:,[1 end]) = true;
    sMerge.imageDensity(:,:,indImage) = ...
        sin(min(bwdist(bound)/sMerge.edgeWidth,1)*pi/2).^2;