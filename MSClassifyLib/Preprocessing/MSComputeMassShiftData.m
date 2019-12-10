function P = MSComputeMassShiftData(msData, varargin)
  % Compute local mass shift profile data
  % P = MSComputeMassShiftData(msData): Compute mass shift profile data
  %   from msData (MSMaldiData object). Result P is an MSFeatureData
  %   object.
  % P = MSComputeMassShiftData(msData, name, value, ...): Specify
  %   additional options as name-value-pairs. Available options are:
  %   - mzRange: Explicit mass range [mzMin, mzMax] (default = [])
  %   - massBins: Number of mass axis intervals (default = 20)
  %   - noiseAdjust: Flag to estimate and adjust for contribution of
  %                  uniform noise (experimental, default = false)
  %   - massScaleDelta: Relative mass defect shift per mass unit
  %                     (default = 4.95e-4, specific to peptides)
  %   - nMoments: Number of moments to compute (experimental, default = 1)
  
  % Check input arguments
  narginchk(1,inf);
  if ~isa(msData, 'MSMaldiData')
    error('Data argument must be an MSMaldiData object');
  end
  
  % Validation functions for optional parameters
  isValidRange = @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 2);
  isScalar = @(x) isscalar(x) && isnumeric(x);
  isScalarInt = @(x) isScalar(x) && x > 0 && mod(x,1) == 0;
  isScalarLogical = @(x) isscalar(x) && islogical(x);
  % Setup parser
  params = inputParser;
  params.addParameter('mzRange', [], isValidRange);
  params.addParameter('massBins', 20, isScalarInt);
  params.addParameter('noiseAdjust', false, isScalarLogical);
  params.addParameter('massScaleDelta', nan, isScalar);
  params.addParameter('nMoments', 1, isScalarInt);
  % Parse input arguments
  params.parse(varargin{:});
  params = params.Results;
  if isempty(params.mzRange)
    params.mzRange = msData.mzVector([1 end]);
  end
  
  % Limit output m/z range to data m/z range
  params.mzRange = [max(params.mzRange(1), msData.mzVector(1)), ...
                    min(params.mzRange(2), msData.mzVector(end))];
  
  % Setup mass defect factor
  if isfinite(params.massScaleDelta)
    mdthRel = params.massScaleDelta;
  else
    mdthRel = 4.95e-4;  % Default value by the peptide mass rule
  end
  mzLambda = 1+mdthRel;

  % Compute distribution profile parameters for all spectra
  [X, mzX] = ComputeMassShift(msData.data, msData.mzVector, params.mzRange, ...
                              params.massBins, mzLambda, params.noiseAdjust, ...
                              params.nMoments);
  X = [real(X), imag(X)];
  valid = reshape(all(isfinite(X), 1), [], 2);
  valid = repmat(all(valid, 2), 1, 2);
  X(:,~valid) = [];
  mzX(~valid(:,1)) = [];
  P = MSFeatureData(X, 'MSComputeMassShiftData', msData);
  P.featureInfo.massShiftProfile = struct;
  P.featureInfo.massShiftProfile.p = params;
  P.featureInfo.massShiftProfile.mzVector = mzX;
end

function [Y, mzY] = ComputeMassShift (X, mzX, mzRange, mzBins, mzLambda, noiseAdjust, nMoments)
  % Compute Y mz bin edges and centers
  indYEdges = interp1(mzX, 1:length(mzX), (0:mzBins)*(mzRange*[-1;1])/mzBins+mzRange(1), 'previous');
  indYEdges(end) = indYEdges(end)+1;
  mzY = (mzX(indYEdges(1:end-1))+mzX(indYEdges(2:end)-1))/2;
  % Compute X bin edges and width
  mzXEdges = [mzX(1:2)*[1.5 -0.5]' (mzX(1:end-1)+mzX(2:end))/2 mzX(end-1:end)*[-0.5 1.5]'];
  dmzX = diff(mzXEdges);
  % Compute X mz weights
  omg = 2*pi/mzLambda;
  c = zeros(nMoments, length(mzX));
  for j = 1:nMoments
    z = exp(1i*j*omg*mzXEdges);
    c(j,:) = -1i/(j*omg)*diff(z)./dmzX;
  end
  if noiseAdjust
    % For noise adjustment, second moments will be computed
    z = exp(2*1i*omg*mzXEdges);
    c2 = -1i/(2*omg)*diff(z)./dmzX;    
  end
  % Initialize Y
  Y = nan(size(X,1), mzBins, nMoments);
  % Iterate over Y bins
  for k = 1:mzBins
    indX = indYEdges(k):indYEdges(k+1)-1;
    Xnorm = sum(X(:,indX),2);
    XnormPos = Xnorm > 0;
    for j = 1:nMoments
      Y(XnormPos,k,j) = X(XnormPos,indX)*c(j,indX).'./Xnorm(XnormPos);
    end
    if noiseAdjust
      % Compute second moment
      Y2 = X(XnormPos,indX)*c2(indX).'./Xnorm(XnormPos);
      % Correction factor to first moment accounting for noise effect
      r = abs(Y(XnormPos,k,1)).^(-4/3) .* abs(Y2).^(1/3);
      % Calculation unstable where first or second moment too small
      rInv = (abs(Y2) < 1e-3) | (abs(Y(XnormPos,k,1)) < 1e-3);
      r(rInv) = 1;
      Y(XnormPos,k,1) = Y(XnormPos,k,1).*r;
    end
    Y(~XnormPos,k,:) = 0;
  end
  Y = reshape(Y, size(Y,1), []);
end
