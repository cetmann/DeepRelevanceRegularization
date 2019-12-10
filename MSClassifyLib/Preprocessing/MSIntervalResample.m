function [msdOut, maskOut] = MSIntervalResample (msdIn, varargin)
  % Resample MALDI data on a specified list of intervals
  %
  % msdOut = MSIntervalResample (msdIn, mzIntervals): Resample msdIn
  %   (MSMaldiData object) on intervals mzIntervals, given as a 2-by-n
  %   matrix of lower and upper interval bounds. The m/z vector of the
  %   resulting msdOut is set to the interval centers.
  %
  % msdOut = MSIntervalResample (msdIn, mzIntervals, mzVector): Explicitly
  %   specify m/z vector of msdOut.
  %
  % msdOut = MSIntervalResample (msdIn, itemMask, mzIntervals, mzVector):
  %   Perform resampling separately for all item subsets specified by
  %   itemMask (MSLabelData, logical vector/matrix, numerical vector). In
  %   case item subsets don't cover all items, msdOut is reduced. Argument
  %   mzIntervals specifies intervals as a cell array of 2-by-n matrices
  %   (see above), where each cell corresponds to one item subset.
  
  % [msdOut, maskOut] = MSIntervalResample(msdIn,mzIntervals): Returns
  % additionally a mask indicating which of the spectra in msdIn are in the
  % msdOut

  % Check input arguments
  narginchk(2,4);
  % Helper functions for argument checking
  isIntervalMatrix = @(x) ~isempty(x) && isnumeric(x) && ismatrix(x) && size(x,1) == 2;
  % Handle different argument list cases
  if nargin == 2
    % Arguments: msdIn, mzIntervals
    itemMask = [];
    mzIntervals = varargin{1};
    mzVector = [];
  elseif nargin == 3
    % Arguments: msdIn, mzIntervals, mzVector
    itemMask = [];
    mzIntervals = varargin{1};
    mzVector = varargin{2};
  else
    % Arguments: msdIn, itemMask, mzIntervals, mzVector
    itemMask = varargin{1};
    mzIntervals = varargin{2};
    mzVector = varargin{3};
  end
  % Convert single intervals matrix to a (1,1) cell array
  if ~iscell(mzIntervals)
    mzIntervals = {mzIntervals};
  end
  
  % Check data argument
  if ~isa(msdIn, 'MSMaldiData')
    error('Data argument must be an MSMaldiData object');
  end
  
  % Check item mask argument
  itemMask = MSCheckItemMask(msdIn, itemMask);
  % Item subsets must not intersect
  itemUnion = sum(itemMask,2);
  if any(itemUnion > 1)
    error('Item subsets must be mutually disjoint');
  else
    % Convert to logical
    itemUnion = itemUnion > 0;
  end
  outNumItems = sum(itemUnion);
  % If item subsets do not cover all data items, output data is reduced
  itemSubsets = find(any(itemMask,1));
  if isempty(itemSubsets)
    error('All item subsets are empty');
  end
  
  % Check intervals argument
  if size(itemMask,2) ~= numel(mzIntervals)
    error('Invalid number of m/z interval matrices specified');
  end
  numIntervals = [];
  for k = itemSubsets
    if ~isIntervalMatrix(mzIntervals{k})
      error('Intervals must be specified as 2-by-n matrix');
    elseif ~all(diff(mzIntervals{k}) > 0)
      error('Intervals must be non-empty (lower < upper bound)');
    elseif min(mzIntervals{k}(:)) < msdIn.mzVector(1) || ...
           max(mzIntervals{k}(:)) > msdIn.mzVector(end)
      error('Intervals must be within data m/z range');
    else
      if isempty(numIntervals)
        numIntervals = size(mzIntervals{k},2);
      elseif size(mzIntervals{k},2) ~= numIntervals
        error('Interval matrices must specify identical number of intervals');
      end
    end
  end
  
  % Check m/z-vector argument
  if ~isempty(mzVector) && ~(isnumeric(mzVector) && isvector(mzVector) && ...
      length(mzVector) == numIntervals)
    error('mzVector must be a numeric vector matching number of intervals');
  elseif isempty(mzVector) && size(itemMask,2) > 1
    error('If itemMask is given, mzVector must be explicitly specified');
  end
    
  % Single item case: By default, new m/z vector is set to interval centers
  if isempty(mzVector)
    mzVector = mean(mzIntervals{1},1);
  end
  % Compute output data
  D = zeros(outNumItems, numIntervals, 'like', msdIn.data);
  % Iterate over non-empty item subsets
  for k = itemSubsets
    D(itemMask(itemUnion,k),:) = ...
      intervalResample(msdIn.data, itemMask(:,k), msdIn.mzVector, mzIntervals{k});
  end
  msdOut = MSMaldiData(D, mzVector, 'normalization', false);
  % Set output data position grid
  if outNumItems == msdIn.numItems
    % No item reduction necessary
    msdOut.setPositions(msdIn.positions);
  elseif ~isempty(msdIn.positions)
    % Reduce position grid
    msdOut.setPositions(msdIn.positions.reduce(itemUnion));
  end
  maskOut = itemUnion;
end


function Y = intervalResample (X, indX, mz, mzI)
  % Sample rows of X indexed by indX with m/z vector mz over intervals
  % mzI (2-by-n matrix), performing linear interpolation at interval bounds.
  % Output Y is the resampled data matrix.

  numI = size(mzI,2);
  Y = zeros(sum(indX), numI, 'like', X);

  % Iterate over intervals
  for k = 1:numI
    % Find part of m/z vector containing the interval, including directly
    % neighboring m/z values left and right of the interval
    leftI = find(mz > mzI(1,k),1);
    rightI = find(mz(leftI-1:end) < mzI(2,k),1,'last');
    indPart = leftI-1:leftI + rightI -1;
    mzPart = mz(indPart);
    w = ones(rightI + 1, 1, 'like', X); % Weighting vector, resampling (~integration)
    % is a weighted sum of corresponding neighborhood masses, those weights
    % will form this w (changing for each mass interval).

    % Considering linear interpolation, resulting in reweighted boundary
    % weights for w.
    mzD1 = (mzPart(2)-mzI(1,k))/(mzPart(2)-mzPart(1));
    mzDe = (mzI(2,k)-mzPart(end-1))/(mzPart(end)-mzPart(end-1));
    w(2) = 0.5+min(0.5,mzD1);
    w(1) = max(0,mzD1-0.5);
    w(end-1) = 0.5+min(0.5,mzDe);
    w(end) = max(0,mzDe-0.5);

    Y(:,k) = X(indX,indPart)*w;
  end
end

