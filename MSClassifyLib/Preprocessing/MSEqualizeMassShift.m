function [msdOut, prfOut] = MSEqualizeMassShift (msdIn, varargin)
  % Equalize mass shifts using mass shift profile normalization
  %
  % msdOut = MSEqualizeMassShift(msdIn): Equalize mass shifts to common
  %   mean mass shift profile
  % msdOut = MSEqualizeMassShift(msdIn, name, value, ...): Specify
  %   additional options for mass shift profile computation as name-value 
  %   pairs (see MSComputeMassShiftData)
  % msdOut = MSEqualizeMassShift(msdIn, prfIn): Equalize mass shifts to
  %   specified reference mass shift profile prfIn (as computed by
  %   MSComputeMassShiftData).
  % msdOut = MSEqualizeMassShift(msdIn, 0): Equalize mass shifts to an
  %   ideal, constant zero reference mass shift profile.
  % msdOut = MSEqualizeMassShift(__, 'unwrap', unwrapFlag): Specify whether
  %   phase unwrapping is applied (default = true).
  % msdOut = MSEqualizeMassShift(__, 'interp', opts): Perform filtering and
  %   interpolation of mass shift profiles before equalization
  %   (experimental, see MSMassShiftInterpolate). Set opts = true to use
  %   default interpolation options, false for no interpolation (default). 
  %   Set opts = {name, value} to pass specific interpolation options.
  % [msdOut, prfOut] = MSEqualizeMassShift(__): Return input data profile
  %   as second argument
  
  % Check input arguments
  narginchk(1,inf);
  if ~isa(msdIn, 'MSMaldiData')
    error('Data argument must be an MSMaldiData object');
  end
  
  % Check whether second argument specifies reference profile (either an
  % MSFeatureData object, the constant zero, or an empty matrix)
  prfIn = [];
  if nargin >= 2
    x = varargin{1};
    if isa(x, 'MSFeatureData')
      if ~isfield(x.featureInfo, 'massShiftProfile')
        error('Invalid input reference mass shift profile');
      end
      xOk = true;
    else
      xOk = isempty(x) || (isnumeric(x) && isscalar(x) && x == 0);
    end
    if xOk
      prfIn = x;
      varargin(1) = [];
    end
  end
  
  % Parse optional name-value parameters
  params = inputParser;
  isLogicalScalar = @(x) isscalar(x) && islogical(x);
  isValidInterp = @(x) isLogicalScalar(x) || iscell(x);
  params.addParameter('unwrap', true, isLogicalScalar);
  params.addParameter('interp', false, isValidInterp);
  params.KeepUnmatched = true;
  params.parse(varargin{:});
  unwrapFlag = params.Results.unwrap;
  interp = params.Results.interp;
  if islogical(interp) && interp
    % If interp is true, use default options
    interp = {};
  end
  
  % Compute mass shift profile of input data
  if isa(prfIn, 'MSFeatureData')
    % Compute profile according to reference profile parameters
    prfParams = {'mzRange', prfIn.featureInfo.massShiftProfile.p.mzRange, ...
                 'massBins', prfIn.featureInfo.massShiftProfile.p.massBins, ...
                 'massScaleDelta', prfIn.featureInfo.massShiftProfile.p.massScaleDelta};
  else
    % No reference profile specified, use specified options
    prfParams = {};
    f = fieldnames(params.Unmatched);
    for k = 1:length(f)
      prfParams = [prfParams, {f{k}, params.Unmatched.(f{k})}];
    end
  end
  prfOut = MSComputeMassShiftData(msdIn, prfParams{:});
  if iscell(interp)
    prfOut = MSMassShiftInterpolate(prfOut, interp{:});
  end

  % If no reference profile specified, use profile of input data
  if isempty(prfIn)
    prfIn = prfOut;
  elseif prfIn == 0
    % If reference profile is specified as 0, create zero shift profile
    m = prfOut.dataLength/2;
    prfIn = MSFeatureData([ones(1,m) zeros(1,m)], 'zero');
    prfIn.featureInfo = prfOut.featureInfo;
  end
  
  % Compute average of reference mass shift profile
  Xm = mean(prfIn.data,1);

  % Transform spectral data
  Y = MassShiftMap(prfOut, Xm, msdIn.data, msdIn.mzVector, unwrapFlag);
  msdOut = MSMaldiData(Y, msdIn.mzVector, 'normalization', false);
  msdOut.setPositions(msdIn.positions);
end

function Y = MassShiftMap (P, M, X, V, W)
  % Apply mass shift map specified by input and output profiles
  % Y = MassShiftMap(P, M, X, V, W)
  %   P: Mass shift profile object corresponding to source data
  %   M: Target mass shift function
  %   X: Source data
  %   V: m/z vector
  %   W: Phase unwrapping flag
  %   Y: Transformed output data
  
  if W
    argFn = @(x,y) unwrap(angle(x+1i*y), [], 2)/(2*pi);
  else
    argFn = @(x,y) mod(angle(x+1i*y)/(2*pi)+0.5, 1)-0.5;
  end
  mzX = P.featureInfo.massShiftProfile.mzVector;
  lx = length(mzX);
  dmzout = argFn(M(1:lx), M(lx+1:end));
  Y = zeros(P.numItems, size(X,2), 'like', X);
  progress = MSProgress('MassShiftMap', P.numItems);
  for k = 1:P.numItems
    progress.update(k);
    dmzin = argFn(P.data(k,1:lx), P.data(k,lx+1:end));
    mzin = interp1(mzX, mzX+dmzin-dmzout, V, 'linear', 'extrap');
    F = griddedInterpolant(V, X(k,:), 'linear', 'nearest');
    Y(k,:) = F(mzin);
  end
  progress.close();
end
