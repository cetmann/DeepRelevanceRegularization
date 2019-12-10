function varargout = MSDataUnion (varargin)
  % Combine multiple data objects into one single object
  % Y = MSDataUnion(X): Combine multiple data objects specified in 2D cell
  %   array X into one single, large output data object Y. If all input
  %   objects include a valid MSPositionGrid object, their spatial
  %   positioning within the output object Y is determined by the layout of
  %   the input cell matrix. If some, but not all input objects include
  %   positions, this information is ignored and a warning is produced.
  %   All non-empty cell entries in X must be of the same object type and
  %   data length. Supported object types include MSData, MSMaldiData, 
  %   MSLabelData, MSFeatureData, and MSDataPartition.
  %
  % [Y1,Y2,...] = MSDataUnion(X1,X2,...): Combine mutltiple sets X1,X2,...
  %   of input objects into output objects Y1,Y2,..., using the same
  %   spatial positioning scheme for all sets. The spatial layout is
  %   determined by X1, position information in all other input sets is
  %   ignored.
  %   Input arguments X1,X2,... must be of same size and filled with data
  %   objects of consistent size at same locations.
  %
  % [__,L,I] = MSDataUnion(__): In addition to the output data objects Yk,
  %   an MSLabelData object L and a cell array I of index vectors may be
  %   returned. Label object L represents the information, which part of
  %   the output data belongs to which input data object. Index vectors I
  %   contain the output object item indices corresponding to each of the
  %   input objects.
  %
  % TODO: Assertions in some loc_ functions
  %       Key-value option to suppress normalization in loc_MaldiUnion
  
  % Check input arguments
  narginchk(1,inf);
  [R,S,T,M] = loc_CellArraysConsistent(varargin{:});
  if ~isempty(R)
    error(['Argument error: ' R]);
  end
  supportedDataTypes = {'MSData', 'MSFeatureData', 'MSMaldiData', 'MSLabelData', 'MSScoreData', 'MSDataPartition'};
  if ~all(ismember(T, supportedDataTypes))
    error('Input arguments contain unsupported data types');
  end
  
  % Number of non-empty input objects
  nM = sum(M(:));
  
  % Check output arguments
  nargoutchk(nargin, nargin+2);
  needOutL = nargout >= nargin+1; % Output argument L requested
  needOutI = nargout >= nargin+2; % Output argument I requested
  
  % Obtain position information from first argument
  P = cell(S);
  % Specify access to input object position information
  if strcmp(T{1}, 'MSDataPartition')
    posFn = @(x) x.classes.positions;
  else
    posFn = @(x) x.positions;
  end
  % Extract position information from input objects in first object set
  P(M) = cellfun(posFn, varargin{1}(M), 'UniformOutput', false);
  
  % Check whether all position objects are present
  nP = sum(~cellfun(@isempty, P(M)));
  if nP < nM
    if nP > 0
      warning('Positions only partly defined, will be ignored');
    end
    P = [];
  end
  
  % Union of position grid objects
  if isempty(P)
    % Positions are not defined, generate item index offsets
    I = loc_ItemIndexOffsets(varargin{1}, M);
    Q = [];
  else
    % Positions are defined, generate union of positions
    [Q,I] = loc_PositionGridUnion(P,M);
  end
  
  % Union of input data object sets, one by one
  for k = 1:nargin
    switch T{k}
      case {'MSData', 'MSFeatureData'}
        varargout{k} = loc_DataUnion(varargin{k}, Q, I, M);
      case 'MSMaldiData'
        varargout{k} = loc_MaldiUnion(varargin{k}, Q, I, M);
      case 'MSLabelData'
        varargout{k} = loc_LabelUnion(varargin{k}, Q, I, M);
      case 'MSScoreData'
        varargout{k} = loc_ScoreUnion(varargin{k}, Q, I, M);
      case 'MSDataPartition'
        varargout{k} = loc_PartitionUnion(varargin{k}, Q, I, M);
      otherwise
        error('Unsupported data type');
    end
  end
  
  % Optional output arguments
  if needOutL
    varargout{nargin+1} = loc_CreateUnionLabel(Q,I,M);
  end
  if needOutI
    varargout{nargin+2} = I;
  end
end


function Y = cellfun_nu (varargin) 
  Y = cellfun(varargin{:}, 'UniformOutput', false);
end


function [R,S,T,M] = loc_CellArraysConsistent (varargin)
  % Check whether input arguments are consistent 2D cell arrays
  % R: Empty if consistent, error message otherwise
  % S: Common size of input cell arrays
  % T: Common data type of non-empty cell entries
  % M: Common mask of non-empty entries
 
  R = [];
  S = [];
  T = [];
  M = [];
  
  % Check whether all arguments are 2D cell arrays
  if ~all(cellfun(@(x) iscell(x) && ismatrix(x), varargin))
    R = 'Not all arguments are 2D cell arrays';
  
  else
    % Check whether all arguments are of same size
    S = cellfun_nu(@size, varargin);
    S = cat(1, S{:});
    if all(min(S,[],1) == max(S,[],1))
      S = S(1,:);
    else
      R = 'Arguments have different sizes';
    end
  end
  
  if isempty(R)
    % Check whether all arguments have non-empty entries at same locations
    M = cellfun_nu(@(x) cellfun(@isempty, x), varargin);
    M = cat(3, M{:});
    if all(all(all(M,3) == any(M,3)))
      M = ~M(:,:,1);
      if ~any(M(:))
        R = 'Arguments must not all be empty';
      end
    else
      R = 'Arguments must have non-empty entries at same locations';
    end
  end
  
  if isempty(R)
    % Check whether non-empty entries within each argument are of same type
    T = cellfun_nu(@(x) unique(cellfun_nu(@class, x(M))), varargin);
    if all(cellfun(@length, T) == 1)
      T = cellfun_nu(@(x) x{1}, T);
    else
      R = 'Non-empty entries within each argument must be of same type';
    end
  end
  
  if isempty(R)
    % Check whether non-empty entries are all scalar
    if ~all(cellfun(@(x) all(cellfun(@isscalar, x(M))), varargin))
      R = 'Non-empty entries must all be scalar';
    end
  end
end


function I = loc_ItemIndexOffsets (P, M)
  % Compute item index offsets
  numItems = zeros(size(P));
  numItems(M) = cellfun(@(p) p.numItems, P(M));
  I = reshape([0 reshape(cumsum(numItems(1:end-1)),1,[])], size(P));
end


function [Q,I] = loc_PositionGridUnion (P, M)
  % Combine multiple position grids into one
  % P: 2D cell array of MSPositionGrid objects (may contain empty cells)
  % M: Logical array indicating non-empty cells
  % Q: Combined MSPositionGrid object
  % I: Matrix of item index offsets corresponding to the input grids P
  
  % Check input arguments
  assert(iscell(P) && ismatrix(P) && islogical(M) && ismatrix(M));
  assert(~isempty(P) && all(size(P) == size(M)));
  assert(all(cellfun(@(x) isa(x, 'MSPositionGrid'), P(M))));
  
  % Compute item index offsets
  I = loc_ItemIndexOffsets(P,M);

  % Get gridSize for all input grids
  pRows = zeros(size(P));
  pCols = zeros(size(P));
  pRows(M) = cellfun(@(p) p.gridSize(1), P(M));
  pCols(M) = cellfun(@(p) p.gridSize(2), P(M));
  % Maximum number of rows/columns for each cell row/column
  nRows = max(pRows, [], 2);
  nCols = max(pCols, [], 1);
  % Cumulative sum of rows/columns, start with an additional zero entry
  cRows = [cumsum(nRows, 'reverse'); 0]; % Assign rows in reverse order
  cCols = [0, cumsum(nCols)];
  
  % Combined indexGrid matrix
  indexGrid = zeros(cRows(1), cCols(end));
  for k = 1:numel(P)
    if ~isempty(P{k})
      [r,c] = ind2sub(size(P),k);
      indexGrid(cRows(r+1)+(1:pRows(r,c)), cCols(c)+(1:pCols(r,c))) = ...
        (P{r,c}.indexGrid+I(k)).*(P{r,c}.indexGrid > 0);
    end
  end
  
  % Output position grid
  Q = MSPositionGrid(indexGrid);
end


function Y = loc_MatrixUnion (X, I, F, M)
  % Concatenate multiple data matrices into one single matrix
  % X: 2D cell array of data matrices (may contain empty cells), must
  %    all have same data type and number of columns
  % I: Matrix of item index offsets corresponding to the input data X, as
  %    returned from PositionGridUnion
  % F: Cell array of function handles corresponding to the input data X
  %    (may contain empty cells). Functions are expected to accept one
  %    argument and to operate on a numeric matrix w/o changing its size.
  % M: Logical array indicating non-empty cells
  % Y: Combined data matrix of same data type as matrices in X.
  
  % Check input arguments
  assert(iscell(X) && ismatrix(X) && ~isempty(X));
  isSameSizeMatrix = @(x) ismatrix(x) && all(size(x) == size(X));
  isCellArrayOfType = @(x,f) all(cellfun(@(y) f(y), x));
  assert(isSameSizeMatrix(M) && islogical(M) && any(M(:)));
  assert(isCellArrayOfType(X(M), @(x) isnumeric(x) && ismatrix(x)));
  assert(isscalar(unique(cellfun(@(x) size(x,2), X(M)))) && ...
         isscalar(unique(cellfun_nu(@(x) class(x), X(M)))));
  assert(isSameSizeMatrix(I) && isnumeric(I) && all(I(:) >= 0) && all(mod(I(:),1) == 0));
  assert(isSameSizeMatrix(F) && isCellArrayOfType(F(M), @(x) isempty(x) || isa(x, 'function_handle')));
  
  nRows = sum(cellfun(@(x) size(x,1), X(M)));
  nCols = size(X{M(1)}, 2);
  dataType = class(X{M(1)});
  Y = zeros(nRows, nCols, dataType);
  for k = 1:numel(X)
    if ~isempty(X{k})
      if isempty(F{k})
        Y(I(k)+(1:size(X{k},1)),:) = X{k};
      else
        Y(I(k)+(1:size(X{k},1)),:) = F{k}(X{k});
      end
    end
  end
end


function Y = loc_DataUnion (X, P, I, M)
  D = cell(size(X));
  D(M) = cellfun_nu(@(x) x.data, X(M));
  F = cell(size(X));
  Y = MSData(loc_MatrixUnion(D,I,F,M));
  Y.setPositions(P);
end


function Y = loc_MaldiUnion (X, P, I, M)
  assert(all(diff(quantile(cell2mat(cellfun_nu(@(x) x.mzVector, reshape(X(M),[],1))), [0 1])) == 0));
  mzv = X{M(1)}.mzVector;
  D = cell(size(X));
  D(M) = cellfun_nu(@(x) x.data, X(M));
  F = cell(size(X));
  Y = MSMaldiData(loc_MatrixUnion(D,I,F,M), mzv);
  Y.setPositions(P);
end


function Y = loc_LabelUnion (X, P, I, M)
  % Combine multiple label data objects into one single object
  % X: 2D cell array of label data objects (may contain empty cells),
  %    must have same number of label columns
  % P: Combined position grid, as returned from PositionGridUnion
  % I: Matrix of item index offsets corresponding to the input array X,
  %    as returned from PositionGridUnion
  % M: Logical array indicating non-empty cells
  % Y: Combined label data object

  % Combine labels
  labels = cellfun_nu(@(x) x.labels(:), X(M));
  labels = unique(cat(1, labels{:}), 'stable');
  
  % Labels to index map for combined labels
  labelsMap = containers.Map(labels, 1:length(labels));
  
  % Old to new index maps and data functions
  indexMaps = cellfun_nu(@(x) [0; cell2mat(labelsMap.values(x.labels(:)))], X(M));
  indexFn = @(ind) @(x) ind(x+1);
  F = cell(size(X));
  F(M) = cellfun_nu(indexFn, indexMaps);
  
  % Combine data
  D = cell(size(X));
  D(M) = cellfun_nu(@(x) x.data, X(M));
  Y = MSLabelData(labels, loc_MatrixUnion(D,I,F,M));
  Y.setPositions(P);
end


function Y = loc_ScoreUnion (X, P, I, M)
  % Combine multiple score data objects into one single object
  % X: 2D cell array of score data objects (may contain empty cells),
  %    must have same labels and number of predictions
  % P: Combined position grid, as returned from PositionGridUnion
  % I: Matrix of item index offsets corresponding to the input array X,
  %    as returned from PositionGridUnion
  % M: Logical array indicating non-empty cells
  % Y: Combined score data object

  % Check that labels are identical in all input objects
  labels = cellfun_nu(@(x) x.labels(:)', X(M));
  assert(isscalar(unique(cellfun(@length, labels))));
  h = cat(1, labels{:});
  assert(size(unique(cellfun(@length, h), 'rows'), 1) == 1);
  assert(size(unique(cell2mat(h), 'rows'), 1) == 1);
  
  % Combine data
  D = cell(size(X));
  D(M) = cellfun_nu(@(x) x.data, X(M));
  F = cell(size(X));
  Y = MSScoreData(labels{1}, loc_MatrixUnion(D,I,F,M));
  Y.setPositions(P);
end


function Y = loc_PartitionUnion (X, P, I, M)
  % Combine multiple partition objects into one single object
  % X: 2D cell array of partition objects (may contain empty cells)
  % P: Combined position grid, as returned from PositionGridUnion
  % I: Matrix of item index offsets corresponding to the input array X,
  %    as returned from PositionGridUnion
  % M: Logical array indicating non-empty cells
  % Y: Combined partition object

  Y = MSDataPartition;
  MM = true(sum(M(:)),1);
  
  % Collect all specific labels
  labels = cellfun_nu(@(x) x.specificLabels(:), X(M));
  labels = unique(cat(1, labels{:}), 'stable');
  for k = 1:length(labels)
    L = cellfun_nu(@(x) x.(labels{k}), reshape(X(M),[],1));
    Y.(labels{k}) = loc_LabelUnion(L, P, reshape(I(M),[],1), MM);
  end
  
  % Collect all other labels
  labels = cellfun_nu(@(x) x.otherLabels(:), X(M));
  labels = unique(cat(1, labels{:}), 'stable');
  for k = 1:length(labels)
    hasLabel = @(x) ~isempty(x) && isfield(x.other, labels{k});
    Mk = cellfun(hasLabel, X(M));
    if all(Mk)
      L = cellfun_nu(@(x) x.other.(labels{k}), reshape(X(M),[],1));
      Y.other.(labels{k}) = loc_LabelUnion(L, P, reshape(I(M),[],1), MM);
    end
  end
end


function L = loc_CreateUnionLabel (P, I, M)
  % Create extra union label indicating indices of original objects
  
  constLabelFn = @(n,k) MSLabelData({num2str(k)}, ones(n,1));
  I = reshape(I(M),[],1);
  numItems = diff([I; P.numItems]);
  ind = find(M(:));
  X = arrayfun(constLabelFn, numItems, ind, 'UniformOutput', false);
  L = loc_LabelUnion(X, P, I, true(length(ind),1));
end



