function [msdOut, prfOut] = MSEqualizeIntensity (msdIn, varargin)
  % Equalize spectral intensities using intensity profile normalization
  %
  % msdOut = MSEqualizeIntensity(msdIn): Equalize intensities to common
  %   mean intensity profile
  % msdOut = MSEqualizeIntensity(msdIn, name, value, ...): Specify
  %   additional options for intensity profile computation as name-value 
  %   pairs (see MSIntensityProfileMap)
  % msdOut = MSEqualizeIntensity(msdIn, prfIn): Equalize intensities to
  %   specified reference intensityProfile prfIn (as computed by
  %   MSIntensityProfileMap)
  % [msdOut, prfOut] = MSEqualizeIntensity(__): Return input data profile
  %   as second argument
  
  % Check input arguments
  narginchk(1,inf);
  if ~isa(msdIn, 'MSMaldiData')
    error('Data argument must be an MSMaldiData object');
  end
  
  % If second argument is an MSFeatureData object, it's the input profile
  prfIn = [];
  if nargin >= 2 && (isa(varargin{1}, 'MSFeatureData') || isempty(varargin{1}))
    prfIn = varargin{1};
    varargin(1) = [];
    if ~isempty(prfIn)
      if ~isfield(prfIn.featureInfo, 'intensityProfile')
        error('Invalid input reference intensity profile');
      elseif ~isempty(varargin)
        error('Input reference profile specified, no additional options supported');
      end
    end
  end
  
  % Compute intensity profile of input data
  if ~isempty(prfIn)
    % Compute profile according to reference profile parameters
    prfParams = {'mzRanges', prfIn.featureInfo.intensityProfile.mzBounds, ...
                 'qScale', prfIn.featureInfo.intensityProfile.qScale};
  else
    % No reference profile specified, use specified options
    prfParams = varargin;
  end
  prfMap = MSIntensityProfileMap(prfParams{:});
  prfOut = prfMap.map(msdIn);

  % If no reference profile specified, use profile of input data
  if isempty(prfIn)
    prfIn = prfOut;
  end
  
  % Compute average of reference intensity profile, ignoring non-finites
  nMz = prfIn.featureInfo.intensityProfile.nMz;
  Xm = reshape(prfIn.data, prfIn.numItems, [], nMz);
  Xm(isinf(Xm(:))) = nan;
  Xm = mean(Xm, 1, 'omitnan');
  % Apply inverse data transformation
  Am = MSIntensityProfileMap.qMapInv(Xm, mean(prfIn.featureInfo.intensityProfile.qData,1));
  % Compute median of distribution over full m/z axis
  qScale = prfIn.featureInfo.intensityProfile.qScale;
  if isvector(Am)
    Atot = mixQuantiles(qScale, Am(:));
  else
    Atot = mixQuantiles(qScale, squeeze(Am));
  end
  md = interp1(qScale, Atot, 0.5);
  % If median is zero, use smallest non-zero quantile value
  if md <= 0
    md = min(Atot(Atot > 0));
  end
  % Normalize reference distribution to median == 1
  Am = Am/md;

  % Transform intensity data
  msdOut = MSMaldiData(IntensityMap(prfOut, Am, msdIn.data, msdIn.mzVector), ...
                       msdIn.mzVector, 'normalization', false);
  msdOut.setPositions(msdIn.positions);
end

function Y = IntensityMap (P, A, X, V)
  % Apply intensity map specified by quantile functions
  % Y = IntensityMap(P, A, X)
  %   P: Intensity profile object
  %   A: Target quantile function
  %   X: Source data
  %   V: m/z vector
  %   Y: Transformed output data
  
  % Get index intervals for m/z ranges
  ind = MSIntensityProfileMap.getIndIntervals(P.featureInfo.intensityProfile.mzBounds, V);
  ind(1,1) = 1;
  ind(2,end) = size(X,2);
  % Number of m/z ranges
  nMz = P.featureInfo.intensityProfile.nMz;
  % Reconstructed quantile functions for each spectrum and range
  B = MSIntensityProfileMap.qMapInv(reshape(P.data, P.numItems, [], nMz), ...
                                    P.featureInfo.intensityProfile.qData);
  % Initialize linear weights to interpolate between m/z ranges
  W = zeros([nMz, size(X,2)], 'like', X);
  for j = 1:nMz
    F = griddedInterpolant(P.featureInfo.intensityProfile.mzCenters, ...
                           1*((1:nMz) == j), 'linear', 'nearest');
    W(j,:) = F(V);
  end
  % Apply per-spectrum intensity transform
  progress = MSProgress('IntensityMap', size(X,1));
  Y = nan(size(X), 'like', X);
  Yk = zeros(size(W), 'like', X);
  for k = 1:size(X,1)
    progress.update(k);
    for j = 1:nMz
      I = ind(1, max(j-1,1)):ind(2, min(j+1,nMz));
      Yk(j,I) = transform(B(k,:,j), A(:,:,j), X(k,I));
    end
    Y(k,:) = sum(W.*Yk);
  end
  progress.close();
end

function Y = transform (B, A, X)
  % Number of leading zeros in A and B
  nzA = find(A > 0, 1)-1;
  nzB = find(B > 0, 1)-1;
  assert(~isempty(nzA), 'Target quantile function must not be all zero');

  if isempty(nzB)
    % Source quantile function all zero, replace last element by median of
    % non-zero input values
    Xgt0 = X > 0;
    if ~any(Xgt0)
      % Input all zero, copy to output and return
      Y = X;
      return
    end
    B(end) = median(X(Xgt0));
    nzB = length(B)-1;
  end
  
  if nzB >= nzA
    % Number of leading zeros determined by B. 
    % Remove leading zeros and add initial (0,0)
    A = [0 A(nzB+1:end)];
    B = [0 B(nzB+1:end)];
  else
    % Number of leading zeros determined by A, remove redundant leading items
    A = A(nzA:end);
    B = B(nzA:end);
  end
  
  % Remove trailing items equal in B
  nzB = find(diff(B) > 0, 1, 'last');
  assert(~isempty(nzB), 'Source quantile function must not be constant');
  A = A(1:nzB+1);
  B = B(1:nzB+1);
  
  % Perform interpolation
  [B,A] = resolveMapAmbiguity(B, A, X, 1e-10);
  F = griddedInterpolant(B, A, 'pchip', 'nearest');
  Y = F(X);
end
