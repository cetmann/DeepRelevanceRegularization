function imsDen = spatialDenoising(ims,foreground,denMethod,denParam)
% Low-level spatial denoising (working on image data)
%
% Input: 
% ims - 3D-array of nz [nx ny] images
%   optional: ...
% foreground - possibility to specify foreground (logical(foreground) will
%           be used in the code), allows boundary consideration so no
%           boundary 0-values will flow onto the foreground
% denMethod - possibility to specify the denoising method with a string
% denParam - possibility to set the denoising parameter(s) in a vector
%
% Output:
% imsDen - 3D-array of nz [nx ny] denoised images
%

narginchk(1,4);
% Set defaults
if nargin < 4 || isempty(denParam)
    denParam = [];
end
if nargin < 3 || isempty(denMethod)
    denMethod = 'convGaussian';
end
if nargin < 2 || isempty(foreground)
    foreground = true(size(ims,1),size(ims,2));
end

foreground = logical(foreground);

switch denMethod
    case {'convOnes','convGaussian'}
        switch denMethod
            case 'convOnes'
                % Set default parameters
                if isempty(denParam)
                    denParam = 2;
                end
                denKernel1D = ones(2*denParam+1,1,'single');
                % Normalization not necessary due to later
                % boundaryCorrection
                % % denKernel = 1/sum(abs(denKernel(:)));
            case 'convGaussian'
                % Set default parameter
                if isempty(denParam)
                    denParam = 2;
                end
                % Determining sigma dependend on denParam (sigma = 1 for 2
                % neighbor pixels)
                sigma = denParam/2;
                gaussianFun = @(x,sigma)(1/sqrt(2*pi*sigma^2)*...
                    exp(-(x.^2/(2*sigma^2))));
                pInt = (-denParam:denParam);
                denKernel1D = single(gaussianFun(pInt',sigma));
                % Normalization not necessary due to later
                % boundaryCorrection                
        end
        
        % Separating convolution filter: 1) direction 1st dim (columns)
        imsDen = convn(ims,denKernel1D,'same');
        % ... 2) direction 2nd dim (rows)
        imsDen = convn(imsDen,denKernel1D','same');
        
        % Boundary Correction
        boundaryCorrection = 1./conv2(denKernel1D,denKernel1D',...
            single(foreground),'same');
        boundaryCorrection(isinf(boundaryCorrection)) = 1;
        boundaryCorrection = repmat(boundaryCorrection,1,1,size(ims,3));
        imsDen = imsDen.*boundaryCorrection;
        
        % Resetting 0s or Nans from background
        imsDen(~foreground) = ims(~foreground);
        
  case 'shock'
    
    denParamToUse = [11,20,0.25,3]; % default
    if ~isempty(denParam)
      denParamToUse(1:length(denParam)) = denParam;
    end
    denParam = denParamToUse;  
    convKernelSize = denParam(1); % Size of convolution kernel (for robust
                                  % gradient approximation)
    numIt = denParam(2); % Number of iterations
    dt = denParam(3); % Time step size of PDE
    preConvolutionSigma = denParam(4);
      
    boundary = imdilate(foreground,true(3,3)) & ~foreground;
    imsDen = spatialDenoising(ims,foreground,'convGaussian',preConvolutionSigma);

    for k = 1:numIt
      %imsDen = spatialDenoising(imsDen,foreground,'convGaussian',1.01);
      imsDenC = spatialDenoising(imsDen,foreground,'convGaussian',convKernelSize);
      imsDenC(isnan(imsDenC)) = 0;
      dilIm = imdilate(imsDenC,true(3,3));
      imsDenC(boundary) = dilIm(boundary);
      [gx,gy] = gradientFunNew(imsDenC);
      normGrad = sqrt(gx.*gx+gy.*gy);
      s = -sign(del2(imsDenC));
      imsDen = imsDen+dt.*s.*normGrad;
      imsDen(~foreground) = ims(~foreground);
    end
    
  otherwise
    error('Unknown spatial denoising method');
end


end
