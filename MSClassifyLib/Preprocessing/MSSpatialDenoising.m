function MSSpatialDenoising(msdata,denMethod,denParam)
% Perform a denoising on an MSData object with the method 'denMethod' and
% parameter 'denParam'. 

% If only msdata is given as an input the default
% denoising method (convGaussian) and parameter will be chosen.
% Uses the low-level function spatialDenoising.m that works on images.
% The meaning of the denParam entries vary from method to method (compare
% spatialDenoising.m for details)
%
% For convolution type denoising, computations are directly performed on
% the data matrix, avoiding the conversion to images.

narginchk(1,3);
if nargin < 3
    denParam = []; % Default parameters will be used
    if nargin < 2 || isempty(denMethod)
        denMethod = 'convGaussian';
    end
end

progress = MSProgress('MSSpatialDenoising', msdata.dataLength);
if any(strcmp(denMethod, {'convOnes', 'convGaussian'}))
    % For convolution type denoising, compute (sparse) convolution matrix
    kernel = spatialKernel(denMethod, denParam);
    A = spatialKernelToMatrix(kernel, msdata.positions);
    % Iterate over m/z positions and apply convolution column-wise
    for k = 1:msdata.dataLength
        progress.update(k);
        % Sparse matrix multiplication does not support single data type,
        % hence convert to double
        msdata.data(:,k) = single(A*double(msdata.data(:,k)));
    end
    
else
    % Use low-level spatialDenoising function for each m/z image
    for k = 1:msdata.dataLength
        progress.update(k);
        msdata.data(:,k) = msdata.positions.decube(...
            spatialDenoising( ...
            msdata.positions.encube(msdata.data(:,k)),...
            msdata.positions.indexGrid,denMethod,denParam)...
            ); 
    end
end
progress.close();

end

