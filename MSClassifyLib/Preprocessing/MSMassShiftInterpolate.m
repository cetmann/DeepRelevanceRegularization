function [Y, R] = MSMassShiftInterpolate (X, varargin)
  % Find and interpolate mass shift profiles violating model assumptions
  %
  % Y = MSMassShiftInterpolate(X): Find spot profiles in mass shift profile
  %   data X where the correlation between decay and m/z is below a
  %   threshold, and replace these by interpolated profiles.
  % Y = MSMassShiftInterpolate(X, name, value, ...): Specify options
  %   'mzRange': m/z range for computing correlations, default=[1500,inf]
  %   'corrType': Correlation type, 'linear' or 'log', default='log'
  %   'corrThresh': Correlation threshold, default=0.9
  %   'interp': Interpolation space, 'profile' or 'pca', default='profile'
  %             Set to 'none' to skip interpolation and return Y empty.
  % [Y,R] = MSMassShiftInterpolate(__): Return MSData object with
  %   correlation coefficients
  
  % Check input arguments
  narginchk(1,inf);
  if ~(isa(X, 'MSFeatureData') && isfield(X.featureInfo, 'massShiftProfile'))
    error('Argument is not a mass shift data object');
  end
  
  % Parse optional name-value parameters
  params = inputParser;
  isPosRange = @(x) isnumeric(x) && isvector(x) && length(x) == 2 && ...
                    all(x > 0) && (diff(x) >= 0);
  isScalar01 = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1);
  isString = @(x) ischar(x) && isvector(x);
  isValidCorrType = @(x) isString(x) && ismember(lower(x), {'linear', 'log'});
  isValidInterp = @(x) isString(x) && ismember(lower(x), {'profile', 'pca', 'none'});
  params.addParameter('mzRange', [1500 inf], isPosRange);
  params.addParameter('corrType', 'log', isValidCorrType);
  params.addParameter('corrThresh', 0.9, isScalar01);
  params.addParameter('interp', 'profile', isValidInterp);
  params.parse(varargin{:});
  
  % Compute profile model correlation coefficients
  profileFn = @(x) x(:,1:end/2)+1i*x(:,end/2+1:end);
  dataFn = @(x,c,f) f(x(:,c));
  switch lower(params.Results.corrType)
    case 'linear'
      corrFn = @abs;
    case 'log'
      corrFn = @(x) log(abs(x));
  end
  mz = X.featureInfo.massShiftProfile.mzVector;
  mzMask = mz >= params.Results.mzRange(1) & mz <= params.Results.mzRange(2);
  R = MSData(-corr(dataFn(profileFn(X.data), mzMask, corrFn)', mz(mzMask)'));
  R.setPositions(X.positions);
  
  % Skip interpolation if 'interp' is set to 'none'
  if strcmpi(params.Results.interp, 'none')
    Y = [];
  else
    % Find spots with correlation below threshold
    M = R.data < params.Results.corrThresh;

    % Perform spatial interpolation
    Y = X.copy;
    switch params.Results.interp
      case 'profile'
        Y.data(M,:) = InterpolateData(X.data, M, X.positions);
      case 'pca'
        C = MSComponentData(X.data, 'mask', ~M);
        C.data(M,:) = InterpolateData(C.data, M, X.positions);
        Y.data(M,:) = C.componentMap(C.data(M,:));
    end
  end
end

function Y = InterpolateData (X, M, P)
  [x,y] = ind2sub(P.gridSize, P.reverseIndex);
  n = size(X,2);
  Y = zeros(sum(M), n);
  for k = 1:n
    F = scatteredInterpolant(x(~M), y(~M), X(~M,k), 'natural', 'nearest');
    Y(:,k) = F(x(M), y(M));
  end
end

