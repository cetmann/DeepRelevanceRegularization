classdef MSPositionGrid < matlab.mixin.Copyable
  % Store a grid of positions associated with spectra in a mass spec dataset
  %
  % This class represents a rectangular grid describing the geometrical
  % positions of data items (spectra, feature vectors, ...) in a mass spec
  % dataset. Non-zero entries in the position grid indicate locations
  % associated with an item, identified by its index within the dataset.
  % The maximum item index is stored as numItems. For consistency reasons,
  % each item index in 1..numItems must occur exactly once in the grid.
  %
  % Properties (read-only):
  %   indexGrid: Index grid where non-zero entries are item indices
  %   gridSize: Size of the position grid
  %   numItems:  Number of items (i.e. spectra) referenced in the grid
  %   reverseIndex: Reverse index vector mapping item indices to grid
  %
  % Methods:
  %   MSPositionGrid: Constructor
  %   encube: Transform matrix of column vectors into 3D cube
  %   decube: Transform 3D cube into matrix of column vectors
  %   reduce: Reduce position grid to specified item subset
  %   compact: Remove stretches of empty space in position grid
  %   assert: Check object consistency, abort if inconsistent
  %
  % MSPositionGrid uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = private)
    indexGrid    % Index grid where non-zero entries are item indices
    numItems     % Number of items (i.e. spectra) referenced in the grid
    
    % Reverse index vector mapping item indices to grid indices 
    % indexGrid(reverseIndex(I)) = I
    %   for I in 1..numItems with reverseIndex(I) ~= nan
    reverseIndex
  end
  
  properties (Dependent)
    gridSize % Size of the position grid

    % Logical vector of length numItems indicating item indices covered by
    % the indexGrid (itemsCovered(I) is true <=> Ex. K: indexGrid(K) = I)
    itemsCovered 
  end
  
  methods
    function obj = MSPositionGrid (indexGrid)
      % Constructor
      % obj.MSPositionGrid(indexGrid): Create position grid. indexGrid must
      %   be a 2D-array of integers [0..numItems] where each non-zero value
      %   occurs at most once.
      
      if ~(isnumeric(indexGrid) && ismatrix(indexGrid) && ~isempty(indexGrid))
        error('indexGrid must be a non-empty 2D integer array');
      else
        obj.indexGrid = indexGrid;
        obj.updateReverseIndex();
      end
    end
    
    function S = get.gridSize (obj)
      S = size(obj.indexGrid);
    end
    
    function M = get.itemsCovered (obj)
      M = ~isnan(obj.reverseIndex);
    end
    
    function N = getMissingItems (obj, numItems)
      % Return number of missing / uncovered items, where numItems is the
      % total number of items in the dataset.
      N = sum(~obj.itemsCovered)+max(numItems-obj.numItems, 0);
    end
    
    function Y = encube (obj, X, bgValue)
      % Transform matrix of column vectors into 3D cube of grid matrices
      % Y = obj.encube(X): Generate 3D cube of grid matrices where the k-th
      %   level Y(:,:,k) represents the data in the k-th column of X.
      % Y = obj.encube(X, bgValue): Fill background with bgValue (default 0)
      
      if ~(isnumeric(X) && ismatrix(X) && size(X,1) == obj.numItems)
        error(['Argument must be a numeric matrix with numItems (= %d)', ...
               'rows'], obj.numItems);
      end
      if nargin >= 3
        if ~(isnumeric(bgValue) && isscalar(bgValue))
          error('Background value must be a numeric scalar');
        end
      else
        bgValue = 0;
      end
      Y = zeros([numel(obj.indexGrid), size(X,2)], 'like', X);
      if bgValue ~= 0
        Y(:) = bgValue;
      end
      for k = 1:size(X,2)
        Y(obj.reverseIndex,k) = X(:,k);
      end
      Y = reshape(Y, [obj.gridSize, size(X,2)]);
    end

    function X = decube (obj, Y)
        % Reverting encube: transform encubed 3D array back into a matrix
        % form with rows sorted corresponding to the indices in the 
        % obj.indexGrid
        
        if ~(isnumeric(Y) && size(Y,1) == size(obj.indexGrid,1) && ...
             size(Y,2) == size(obj.indexGrid,2))
            error(['Argument must be a numeric array with the first ', ...
               '2 dimensions sizes matching the indexGrid ([%d %d])'], ...
               size(obj.indexGrid,1),size(obj.indexGrid,2));
        end        
        [nx,ny,nz] = size(Y);
        Y = reshape(Y,nx*ny,nz);
        X = Y(obj.reverseIndex,:);   
    end

    function [P, offset] = reduce (obj, itemMask)
        % Reduce the MSPositionGrid object to a version constrained to the
        % logical item Mask vector
        %
        % Note: this method does not allow to modify the object but
        % only generates a new object. 
        % Otherwise it would be very easy to destroy the structure of 
        % MSData objects by typing msData.positions.reduce(itemMask).
        % If an entire data object has to be reduced use the MSData.reduce 
        % method instead.
        %
        % The indices of are monotonically rewritten
        % such that...
        % ... data objects generated as a reduction onto itemMask (eg.
        %   the pseudoImages in LinearBasisMaps determined only on a 
        %   training part of the data) can be visualized directly when they
        %   get assigned a reduced position grid version that can be
        %   generated with this method.
        
        % Check input argument
        if ~((isnumeric(itemMask) || islogical(itemMask)) && ...
             isvector(itemMask) && length(itemMask) == obj.numItems)
          error(['itemMask must be a numeric or logical vector of length ' ...
                 'obj.numItems (=%d)'], obj.numItems);
        elseif ~any(itemMask)
          error('itemMask is empty, no items selected');
        end
        
        % Generate map of old to new grid indices
        map = zeros(1, obj.numItems+1);
        selectedItems = find(itemMask);
        map(selectedItems+1) = 1:length(selectedItems);
        % Generate position grid with reduced item indices
        grid = map(obj.indexGrid+1);
        % Crop grid to bounding box
        [grid, offset] = crop2d(grid);
        % Create a new MSPositionGrid object
        P = MSPositionGrid(grid);   
    end
    
    function compact (obj, N)
      % Remove stretches of empty space in position grid
      % obj.compact(): Replace stretches of two or more empty rows and
      %   columns in position grid by single empty rows / columns
      % obj.compact(N): Specify minimum length of empty row / column
      %   stretches (default: N=1)
      % obj.compact([N1, N2]): Specify minimum stretch length for rows (N1)
      %   and columns (N2) separately
      % obj.compact('pack'): Remove all empty rows and columns
      % obj.compact(M): Remove all empty rows and columns that are not
      %   covered by mask M (logical matrix matching grid size)
      
      % Check arguments
      M = [];
      if nargin < 2 || isempty(N)
        N = 1;
      elseif ischar(N) && isvector(N) && strcmp(N, 'pack')
        N = [];
        M = false(obj.gridSize);
      elseif islogical(N)
        if ~(ismatrix(N) && all(size(N) == obj.gridSize))
          error('Item mask must match grid size');
        end
        M = N;
        N = [];
      elseif ~(isnumeric(N) && isvector(N) && length(N) <= 2 && ...
               all(N >= 1) && all(mod(N,1) == 0))
        error('N must be a positive integer scalar or vector of length 2');
      end
      if isscalar(N)
        N = [N N];
      end

      % Iterate over dimensions
      M0stretch = cell(1,2);
      for d = 1:2
        if isempty(N)
          M0stretch{d} = all(obj.indexGrid == 0 & ~M, d);
        else
          % Require at least N+1 lines (rows/columns)
          if size(obj.indexGrid,d) > N(d)
            % Find all-zero lines
            M0 = all(obj.indexGrid == 0, d);
            % Find starting positions of stretches of N+1 all-zero lines
            M0stretch{d} = [all(hankel(M0(1:N(d)+1), M0(N(d)+1:end)),1) ...
                            false(1,N(d))];
          end
        end
      end
      % Remove stretches
      obj.indexGrid(:,M0stretch{1}) = [];
      obj.indexGrid(M0stretch{2},:) = [];
      % Update reverse index, checking for inconsistencies
      prevCovered = obj.itemsCovered;
      obj.updateReverseIndex();
      newCovered = obj.itemsCovered;
      assert(all(size(prevCovered) == size(newCovered)) && all(prevCovered == newCovered));
    end

    function assert (~)
      % Check object consistency, abort if inconsistent
      
      % Implemented empty, access restrictions enforce consistency
    end
    
  end
  
  methods (Access = protected)
    function updateReverseIndex (obj)
      % Compute reverseIndex and numItems from indexGrid
      
      % Sort indexGrid entries
      s = zeros(numel(obj.indexGrid),2);
      [s(:,1),s(:,2)] = sort(obj.indexGrid(:));
      % Reduce to non-zero entries
      s = s(s(:,1) ~= 0,:);
      % Check that entries are unique
      d = diff(s(:,1));
      if any(d == 0)
        error('Non-zero entries of indexGrid must be unique');
      end
      % Return results
      obj.numItems = s(end,1);
      obj.reverseIndex = nan(obj.numItems,1);
      obj.reverseIndex(s(:,1)) = s(:,2);
    end
  end
  
end


% local functions 

function [Y, offset] = crop2d (X)
  % Crop 2D array to bounding box of non-zero elements
  bbox = nan(2,2);
  for d = 1:2
    nonZero = find(any(X,d));
    if ~isempty(nonZero)
      bbox(d,:) = nonZero([1 end]);
    end
  end
  if any(isnan(bbox))
    Y = [];
    offset = [0 0];
  else
    Y = X(bbox(2,1):bbox(2,2), bbox(1,1):bbox(1,2));
    offset = bbox([2 1],1);
  end
end

