function kernel = spatialKernel (kernelType, kernelParam)
  % Compute convolution kernel for spatial denoising
  %
  % Input: 
  % kernelType - Kernel type (default is 'convGaussian')
  %   'convOnes':     Square, constant kernel
  %   'convGaussian': Square gaussian kernel with sigma = kernelParam/2
  % kernelParam - Kernel radius (default is 2).  Kernel size is 
  %               (2*kernelParam+1) squared.
  %
  % Output:
  % kernel - Kernel matrix
  %

  narginchk(0,2);
  % Set defaults
  if nargin < 2 || isempty(kernelParam)
    kernelParam = 2;
  end
  if nargin < 1 || isempty(kernelType)
    kernelType = 'convGaussian';
  end

  switch kernelType
    case 'convOnes'
      denKernel1D = ones(2*kernelParam+1,1);
    case 'convGaussian'
      % Determining sigma dependend on denParam (sigma = 1 for 2
      % neighbor pixels)
      sigma = kernelParam/2;
      gaussianFun = @(x,sigma)(1/sqrt(2*pi*sigma^2)*exp(-(x.^2/(2*sigma^2))));
      pInt = (-kernelParam:kernelParam);
      denKernel1D = gaussianFun(pInt', sigma);

    otherwise
      error('Unknown kernel type');
  end
  
  % Compute 2D kernel and normalize to sum == 1
  kernel = denKernel1D*denKernel1D';
  kernel = kernel/sum(kernel(:));
end
