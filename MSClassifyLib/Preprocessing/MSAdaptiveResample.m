function [msdOut, maskOut] = MSAdaptiveResample (msdIn, varargin)
  % Perform adaptive resampling based on the data mass shift profile
  %
  % msdOut = MSAdaptiveResample(msdIn): Perform adaptive resampling of
  %   msdIn to m/z intervals derived from integer nominal masses
  %
  % msdOut = MSAdaptiveResample(msdIn, 'equalize'): Use mass shift
  %   equalization to process each spectrum separately
  %
  % msdOut = MSAdaptiveResample(msdIn, itemMask): Perform adaptive
  %   resampling separately for all item subsets specified by itemMask
  %   (MSLabelData, logical vector/matrix, numerical vector). In case item
  %   subsets don't cover all items, msdOut is reduced accordingly.
  %
  % msdOut = MSAdaptiveResample(msdIn, itemMask, profile): Use given mass
  %   shift profile instead of computing it from msdIn. Number of profile
  %   items and item subsets must match. Use empty itemMask to select all
  %   items.
  %
  % msdOut = MSAdaptiveResample(__, name, value, ...): Specify additional
  %   parameters as name value pairs (see MSAdaptiveResampleParams)
  
  % [msdOut, maskOut] = MSAdaptiveResample(msdIn, varargin): Returns
  % additionally which spectra from the msdIn are in msdOut

  % Check input arguments
  narginchk(1,inf);
  equalize = false;
  itemMask = [];
  profile = [];
  if nargin >= 2
    if ~ischar(varargin{1})
      % If second argument is not a string, it's assumed to be itemMask
      itemMask = varargin{1};
      varargin(1) = [];
      if nargin >= 3 && ~ischar(varargin{1})
        % If third argument is not a string, it's assumed to be profile
        profile = varargin{1};
        varargin(1) = [];
      end
    elseif strcmpi(varargin{1}, 'equalize')
      % If second argument is 'equalize', select mass shift equalization
      equalize = true;
      varargin(1) = [];
    end
  end
  if ~isa(msdIn, 'MSMaldiData')
    error('Data argument must be an MSMaldiData object');
  end
  itemMask = MSCheckItemMask(msdIn, itemMask);
  if ~isempty(profile) && profile.numItems ~= size(itemMask,2)
    error('Number of profile items and item subsets must match');
  end
  % Optional parameters
  P = MSAdaptiveResampleParams(varargin{:});
  % mzType 'center' not allowed with multiple item subsets
  if strcmp(P.mzType, 'center') && size(itemMask,2) > 1 || ...
     (~isempty(profile) && profile.numItems > 1)
    error('mzType ''center'' not allowed with multiple item subsets');
  end
  
  % If requested, perform mass shift equalization on input data
  if equalize
    msdIn = MSEqualizeMassShift(msdIn, 'mzRange', P.mzRange, ...
                                'massScaleDelta', P.massScaleDelta, ...
                                'unwrap', P.unwrap);
  end
  
  % Generate mass shift profile
  if isempty(profile)
    profile = MSMassShiftProfile('massScaleDelta', P.massScaleDelta, 'unwrap', P.unwrap);
    profile.computeProfile(msdIn, itemMask);
  end
  
  % m/z range and vector
  if isempty(P.mzRange)
    % Compute from profile, requiring peptide signal in at least half of
    % the profile items
    P.mzRange = profile.signalRange(ceil(profile.numItems/2));
  else
    % Limit m/z range to intersection of specified and original range
    P.mzRange = [max(P.mzRange(1), profile.mzRange(1)), ...
                 min(P.mzRange(2), profile.mzRange(2))];
  end
  % Nominal mass range and vector
  mzRangeNom = profile.exactToNominalMass(P.mzRange);
  mzNom = mzRangeNom(1):(1/P.subsampling):mzRangeNom(2);
  % Exact mass vector
  mzExact = profile.nominalToExactMass(mzNom);
  
  % Apply mass shifts to exact mass vector and compute intervals
  mzCenter = repmat(mzExact', [1 profile.numItems])+profile.massShift(mzExact');
  mzIntervals = cat(3, mzCenter-P.width/2, mzCenter+P.width/2);
  % Make sure m/z intervals lie within data m/z range
  validIntervals = all(mzIntervals(:,:,1) >= msdIn.mzVector(1), 2) & ...
                   all(mzIntervals(:,:,2) <= msdIn.mzVector(end), 2);
  mzIntervals = mzIntervals(validIntervals,:,:);
  % Convert into cell array of interval matrices
  mzIntervals = squeeze(num2cell(permute(mzIntervals, [3 1 2]), [1 2]));
  
  % Apply interval resample
  switch P.mzType
    case 'exact'
      mzVector = mzExact(validIntervals);
    case 'nominal'
      mzVector = mzNom(validIntervals);
    case 'center'
      assert(isvector(mzCenter));
      mzVector = mzCenter(validIntervals)';
    otherwise
      error('Invalid mzType');
  end
  [msdOut, maskOut] = MSIntervalResample(msdIn, itemMask, mzIntervals, mzVector);
end

