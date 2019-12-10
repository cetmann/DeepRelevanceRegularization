function Y = MSMeanAndBounds (X, dim, boundsType)
  % Return mean and error bounds of an input data array
  % Y = MSMeanAndBounds(X): Compute mean and difference to minimum and 
  %   maximum of data values in X along the last dimension of X and
  %   concatenate them in that dimension.
  % Y = MSMeanAndBounds(X, dim): Compute mean and bounds along dimension
  %   given in dim (positive integer scalar) and return results concatenated
  %   along this dimesion.
  % Y = MSMeanAndBounds(X, dim, boundsType): Specify type of bounds to
  %   compute:
  %   - 'range': Minimum and maximum values (default)
  %   - 'std'  : Standard deviation
  %   - Q      : Quantiles Q and 1-Q (Q: real scalar in range 0..1)
  %   Specify [] as dim argument to compute mean and bounds along the last
  %   dimension of X.

  % Check input arguments
  narginchk(1,3)
  if nargin < 2
    dim = [];
  end
  if nargin < 3
    boundsType = 'range';
  end
  if ~(isnumeric(X) && ~isempty(X))
    error('Data argument X must be a non-empty numeric array');
  elseif ~(isnumeric(dim) && (isempty(dim) || (isscalar(dim) && dim > 0)))
    error('Dimension argument must be a positive integer scalar');
  elseif ~(ischar(boundsType) || (isnumeric(boundsType) && ...
           isscalar(boundsType) && boundsType >= 0 && boundsType <= 1))
    error ('Bounds type must be a string value or a real scalar in range 0..1')
  end
  
  % If dim is not specified, operate on the last dimension of X
  if isempty(dim)
    dim = ndims(X);
  end
  
  % Compute mean of X
  M = mean(X, dim);
  
  % Compute bounds according to boundsType
  if strcmp(boundsType, 'std')
    S = std(X, 0, dim);  % Standard deviation, weighted 1/(N-1)
    Y = cat(dim, M, S, S);
  else
    if strcmp(boundsType, 'range')
      % Difference between mean and minimum / maximum
      Y = cat(dim, M, M-min(X,[],dim), max(X,[],dim)-M);
    elseif isnumeric(boundsType)
      Q = quantile(X, sort([boundsType 1-boundsType]), dim);
      Y = cat(dim, M, abs(cat(dim,M,M)-Q));
    else
      error('Invalid bounds type specified');
    end
  end
end

