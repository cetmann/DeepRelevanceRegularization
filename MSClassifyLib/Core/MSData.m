 classdef MSData < matlab.mixin.Copyable
  % Base class for mass spec data
  %
  % This class is a base class for storing mass spec datasets, including
  % original spectral data as well as derived feature data. Data is stored
  % as a 2D-array where rows represent data items (i.e. spectra or
  % feature vectors). Additional, optional properties may describe the
  % mz-values related to the columns of the data matrix, annotation
  % information (labels or regions) and position information (geometrical
  % layout).
  %
  % Properties:
  %   data: Data matrix, items (spectra, feature vectors) stored as rows
  %   dataLength: Length of data items (number of mz-values or features)
  %   numItems: Number of data items
  %   annotations: MSAnnotationSet object (optional)
  %   positions: MSPositionGrid object (optional)
  %
  % Methods:
  %   MSData: Constructor
  %   setAnnotations: Set annotation set object
  %   setPositions: Set position grid object
  %   partition: Create an MSLabelData object representing a data partition
  %   reduce: Reduce data items to specified item subset
  %   meanData: Compute mean data item
  %   show: Show images of data at selected indices
  %   assert: Check object consistency, abort if inconsistent
  %   getChunks: Static utility method to subdivide [1..N] into chunks
  %
  % MSData uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a (shallow) copy.
  
  properties
    % Data matrix, data items stored as rows
    % Rows of data represent data items, such as spectra or feature 
    % vectors. Data matrix is either initialized in constructor or via 
    % assignment.
    % Data is available for public access in order to simplify manipulation
    % of individual items and avoid implicit copying of data. However,
    % changing the matrix size may result in incosistencies with other
    % object properties.
    data = [];
    additionalDataVersions = {};
  end
  properties (SetAccess = protected)
    annotations = []; % MSAnnotationSet object (optional)
    positions = [];   % MSPositionGrid object (optional)
  end
  properties (Dependent)
    dataLength % Length of data items (number of mz-values or features)
    numItems   % Number of data items (spectra)
  end
  
  methods
    function obj = MSData (varargin)
      % Constructor
      % obj = MSData: Create an empty MSData object. Data must be
      %   initialized by set access.
      % obj = MSData(data): Create an MSData object for the given data
      %   matrix. Optional properties (annotations, position grid) are
      %   initialized empty.

      narginchk(0,1);
      if nargin == 1
        if ~(ismatrix(varargin{1}) && isnumeric(varargin{1}) && ~isempty(varargin{1}))
          error('Data matrix must be a non-empty, numeric 2D-array');
        else
          obj.data = varargin{1};
        end
      end
    end
    
    function N = get.dataLength (obj)
      % Get length of data items (= number of columns in data matrix)
      N = size(obj.data,2);
    end
    
    function N = get.numItems (obj)
      % Get number of data items (= number of rows in data matrix)
      N = size(obj.data,1);
    end
    
    function setAnnotations (obj, A)
      % Set annotation set object
      % obj.setAnnotations(A): Set MSAnnotationSet object A. Numbers of data
      %   items must match.
      
      if ~isempty(A) && ~(isa(A, 'MSAnnotationSet') && A.numItems == obj.numItems)
        error(['annotations must be an MSAnnotationSet object with numItems ' ...
               'matching the number of data rows']);
      else
        obj.annotations = A;
      end
    end
    
    function setPositions (obj, P)
      % Set position grid object
      % obj.setPositions(P): Set MSPositionGrid object P. Numbers of data
      %   items must match.
      if ~isempty(P) && ~(isa(P, 'MSPositionGrid') && P.numItems == obj.numItems)
        error(['positions must be an MSPositionGrid object with numItems ' ...
               'matching the numer of data rows']);
      else
        obj.positions = P;
      end
    end
    
    function Y = meanData (obj, varargin)
      % Compute mean data
      % Y = obj.meanData: Compute mean of all data items, result is a row
      %   vector
      % Y = obj.meanData(A): Compute means of data items included in 
      %   specified annotations, given as a vector of index values or as a
      %   regular expression name pattern. Means are computed separately for
      %   each annotation and are returned as rows of the output matrix Y.
      % Y = obj.meanData(L): For MSLabelData object L, compute means of
      %   data items grouped by non-zero label values in L
      % Y = obj.meanData(M): For logical item mask M (vector or matrix), 
      %   compute means of data items corresponding to columns of M
      % Y = obj.meanData(__, options): Specify additional options using
      %   keywords. The following are supported
      %   - 'includenan', 'omitnan' (default): Include or omit NaN values
      %   - 'sum', 'mean' (default): Compute sum or mean of data
      
      narginchk(1,4);
      % Specify optional flag strings (default is first)
      optFlags = {{'omitnan', 'includenan'}; ...
                  {'mean', 'sum'}};
      % Specify corresponding option values
      optValues = {optFlags{1}; ...
                   {@mean, @sum}};
      % Parse optional flag arguments
      [optInds, varargin] = obj.parseFlagOptions(varargin, optFlags);
      % Assign selected (or default) values
      selectedValue = @(k) optValues{k}{max(optInds(k), 1)};
      nanflag = selectedValue(1);
      accuFn = selectedValue(2);
      
      if length(varargin) > 1
        error('Invalid options specified');
      elseif isempty(varargin)
        % No annotations specified, compute mean over all data items
        Y = accuFn(obj.data, 1, nanflag);
      else
        A = varargin{1};
        if isa(A, 'MSLabelData')
          % Label data specified, numItems must match
          if obj.numItems ~= A.numItems
            error('Number of items of label argument does not match');
          end
          M = A.items;
        elseif islogical(A)
          if ~ismatrix(A)
            error('Item mask argument must be a vector or matrix');
          end
          if isvector(A)
            A = A(:);
          end
          if obj.numItems ~= size(A,1)
            error('Size of item mask argument does not match');
          end
          M = A;
        else
          % Annotations specified
          if isempty(obj.annotations)
            error('No annotations specified for this data object');
          else
            M = obj.annotations.items(A);
          end
        end
        Y = zeros(size(M,2), obj.dataLength, 'like', obj.data);
        for k = 1:size(M,2)
          Y(k,:) = accuFn(obj.data(M(:,k),:), 1, nanflag);
        end
      end
    end
    
    function labelData = partition (obj, labels, annotations)
      % Create an MSLabelData object representing a data partition
      % labelData = obj.partition(labels, annotations):
      %   Create an MSLabelData object representing the partition of the
      %   data items into disjoint segments
      %   labels: Cell string vector of names associated with segments
      %   annotations: Cell vector of annotations specifiying disjoint 
      %     segments. Each cell element may be vector of annotation indices
      %     or a regular expression annotation name pattern.
      %
      %   labels and annotations must be cell vectors of same length. If
      %   only one label is defined, both 1-by-1 cell arrays may be
      %   replaced by the respective cell elements.
      %
      %   The returned labelData is an MSLabelData object with label values
      %   in range 0..number of labels. Postive label values correspond to
      %   specified labels, 0 corresponds to the complement of the union of
      %   all labels.
      
      % Assert consistency
      obj.assert;
      % If arguments are not cell arrays, convert to 1-by-1 cell arrays
      if ~iscell(labels)
        labels = {labels};
      end
      if ~iscell(annotations)
        annotations = {annotations};
      end
      % Check arguments
      if ~(isNonEmptyCellVector(labels) && isNonEmptyCellVector(annotations) && ...
           length(labels) == length(annotations))
         error('labels and annotations must be cell array vectors of same length')
      elseif ~iscellstr(labels)
        error('labels must be a cell string array')
      end
      
      % Create label data from annotations
      labelData = MSLabelData(labels, obj.annotations.segmentLabels(annotations), obj);
    end
    
    function varargout = reduce (obj, S, keepNumItemsFlag)
      % Reduce data items to specified item subset
      % obj.reduce(S): Reduce data items to the subset M specified by S.
      %   S may be one of the following:
      %     - a boolean vector of length == numItems; is copied to M
      %     - a regular expression name pattern or an index vector selecting
      %       a set of annotations; M is the union of all annotations
      %     - an MSLabelData object; M is the set of all items with non-zero
      %       label values
      %   If obj includes an annotation set and/or a position grid, these
      %   are reduced consistently.
      % D = obj.reduce(S): Same as above, but create and return a new data
      %   object of the same type as obj. The object obj itself is not
      %   changed.
      % If keepNumItemsFlag is specified as true (default is false), the
      % reduction effectively sets ~S to zero but doesn't change the size
      % of the object, i.e. numItems will remain the same and the position
      % and annotation reduce methods will not be applied. 
      % HINT: This option is useful mainly for label data, and might lead 
      % to undesired results when used carelessly for MSMaldiData objects!
      
      % Assert consistency
      obj.assert;
      % Check number of output and input arguments
      nargoutchk(0,1);
      narginchk(2,3);
      % Convert item selector to boolean vector
      itemMask = obj.resolveItemSelector(S);
      
      if nargin < 3
        keepNumItemsFlag = false;
      end

      % Selected item subset empty?
      if ~any(itemMask)
        if nargout >= 1
          D = [];
        else
          error('Cannot reduce to empty item subset');
        end
        return;
      end
      
      % Create reduced annotation set if present
      if isempty(obj.annotations)
        A = [];
      elseif keepNumItemsFlag
        A = obj.annotations;
      else
        A = obj.annotations.reduce(itemMask);
      end
      
      % Create reduced position grid if present
      if isempty(obj.positions)
        P = [];
      elseif keepNumItemsFlag
        P = obj.positions;
      else        
        P = obj.positions.reduce(itemMask);
      end
      
      % Modify object in place or create new object?
      
      if ~keepNumItemsFlag
        if nargout >= 1
          % Create a new object of the same type as obj with the reduced data
          D = obj.newFromTemplate(obj.data(itemMask,:));
        else
          % Modify data
          obj.data = obj.data(itemMask,:);
          D = obj;
        end
      else
        dataMatrixSet0 = obj.data;
        dataMatrixSet0(~itemMask,:) = 0;
        if nargout >= 1
          D = obj.newFromTemplate(dataMatrixSet0);
        else
          obj.data = dataMatrixSet0;
          D = obj;
        end
      end
      
      % Associate new annotations and position grid
      if ~isempty(A)
        D.setAnnotations(A);
      end
      if ~isempty(P)
        D.setPositions(P);
      end
      
      % Return new object D if requested
      if nargout >= 1
        varargout{1} = D;
      end
    end
    
    function show (obj, I, varargin)
      % Show images of data at selected indices
      % obj.show(I): Create figure with image plots of data at indices
      %   specified in I
      % obj.show: Plot first 6 images
      % obj.show(I, Name, Value, ...): Specify additional parameters as
      %   name-value pairs. Supported parameters are:
      %   target: 'figure'  - Create new figure (default)
      %           'current' - Plot in current figure/axes
      %           'clear'   - Clear current figure before plotting
      
      % Assert consistency
      obj.assert;
      % Check arguments
      if isempty(obj.positions)
        error('No positions specified for this data object, cannot show')
      end
      if nargin < 2 || isempty(I)
        % No data indices specified, show the first six
        I = 1:min(obj.dataLength, 6);
      end
      if ~(isnumeric(I) && isvector(I) && ~isempty(I) && min(I) >= 1 && ...
           max(I) <= obj.dataLength)
         error(['Index argument must be an integer vector with values in' ...
                'range 1..dataLength (= %d)'], obj.dataLength)
      end
      
      % Parse parameters
      params = inputParser;
      isValidTarget = @(x) ismember(x, {'figure', 'current', 'clear'});
      params.addOptional('target', 'figure', isValidTarget);
      params.parse(varargin{:});
      if strcmp(params.Results.target, 'current') && length(I) > 1
        error('Cannot plot multiple images with target=current');
      end
      
      % Setup figure
      tiles = MSFitTiling(length(I), 4/3);
      switch params.Results.target
        case 'figure'
          figure
        case 'clear'
          clf
      end
      % Plot images
      for k = 1:length(I)
        if ~strcmp(params.Results.target, 'current')
          subplot(tiles(1), tiles(2), k);
        end
        obj.showImage(obj.data(:,I(k)));
      end
    end
      
    function assert (obj)
      % Check object consistency, abort if inconsistent
      
      % Check validity of data member
      assert(ismatrix(obj.data) && isnumeric(obj.data), ...
             'Assertion failed: MSData.data is not a numeric matrix');
      % Check consistency with annotations member
      if ~isempty(obj.annotations)
        obj.annotations.assert();
        assert(size(obj.data,1) == obj.annotations.numItems, ...
               ['Assertion failed: MSData.data row count (%d) ', ...
                'does not match annotations.numItems (%d)'], ...
               size(obj.data,1), obj.annotations.numItems);
      end
      if ~isempty(obj.positions)
        obj.positions.assert();
        assert(size(obj.data,1) == obj.positions.numItems, ...
               ['Assertion failed: MSData.data row count (%d) ', ...
                'does not match positions.numItems (%d)'], ...
               size(obj.data,1), obj.positions.numItems);
      end
    end
  
  end
  
  methods (Access = protected)
    function C = getColorMap (~)
      % Return the color map to be used by the show method
      % May be overloaded in derived classes
      
      % Use virids with background color (value == 0) set to black
      C = [zeros(1,3); viridis(255)];
    end
    
    function showImage (obj, X)
      % Show vector X as image
      foreground = logical(obj.positions.indexGrid);
      imagesc(flipud((obj.positions.encube(X))), ...
              'alphadata', flipud(foreground), 'clipping','off');
      axis image; axis off;
      colormap(gca, obj.getColorMap);
    end
    
    function M = resolveItemSelector (obj, S)
      % Return boolean vector converted from item selector argument S.
      % S may be one of the following:
      %   - a boolean vector of length == numItems; is copied to M
      %   - a regular expression name pattern or an index vector selecting
      %     a set of annotations; M is the union of all annotations
      %   - an MSLabelData object; M is the set of all items with non-zero
      %     label values
      if islogical(S)
        if ~(isvector(S) && length(S) == obj.numItems)
          error(['Boolean argument I must be a vector of length matching ' ...
                 'numItems (= %d)'], obj.numItems);
        end
        M = S;
      elseif isa(S, 'MSLabelData')
        if ~(S.numItems == obj.numItems)
          error('MSLabelData argument I must match numItems (= %d)', obj.numItems);
        end
        M = any(S.data > 0, 2);
      elseif isempty(obj.annotations)
        error('No annotations specified for this data object');
      else
        M = any(obj.annotations.items(S), 2);
      end
    end
    
    function D = newFromTemplate (~, data)
      % Create a new MSData object with specified data
      % This method must be overloaded by derived classes to create objects
      % of the derived class and possibly to copy properties from obj to
      % the new object.
      D = MSData(data);
    end
  end
  
  methods (Static)
    function Y = getChunks (N, Kmax)
      % Subdivide integer range 1:N into a list of K intervals, where K is
      % chosen as the minimum of sqrt(N) and Kmax (default = 10).
      % Result Y is a K-by-2 matrix of integers
      
      narginchk(1,2)
      if nargin < 2
        Kmax = [];
      end
      if isempty(Kmax)
        Kmax = 10;
      end
      
      K = min(Kmax, ceil(sqrt(N)));
      L = ceil(N/K);
      Y = [(0:K-1)'*L+1, min((1:K)'*L, N)];
    end
  end
  
  methods (Static, Access = public)
    function [opts, args] = parseFlagOptions (args, optStrings)
      % Search argument list args for optional flag strings specified in
      % optStrings. Return indices of matched flags in opts, remaining
      % arguments in args.
      
      % For each list of flags in optStrings ...
      opts = zeros(length(optStrings),1);
      for k = 1:length(optStrings)
        % Find matching string arguments in args
        iStrArgs = find(cellfun(@ischar, args));
        iOptArgs = iStrArgs(ismember(args(iStrArgs), optStrings{k}));
        if length(iOptArgs) > 1
          error('Only one of the options %s may be specified', strjoin(optStrings{k}, ', '));
        elseif length(iOptArgs) == 1
          opts(k) = find(strcmp(args{iOptArgs}, optStrings{k}));
          args(iOptArgs) = [];
        end
      end
    end
  end
    
end


function result = isNonEmptyCellVector (X)
  result = iscell(X) && isvector(X) && ~isempty(X);
end

