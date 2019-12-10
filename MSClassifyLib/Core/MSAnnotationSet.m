classdef MSAnnotationSet < matlab.mixin.Copyable
  % Store and access annotations for mass spec data
  %
  % This class represents a set of annotations for mass spec datasets 
  % consisting of an ordered collection of data items (spectra, feature 
  % vectors, etc.). 
  % Each annotation is a subset of items associated with a name. The
  % subsets are represented as logical vectors (masks) of length equaling
  % the number of data items.
  % Annotations may be addressed via their index or via their name. Regular
  % expressions may be used to select multiple annotations matching a name 
  % pattern.
  %
  % Properties (read-only):
  %   numItems: Total number of items in the dataset
  %   annotations: List of annotations with name and mask vector
  %   numAnnotations: Number of annotations
  %   attributes: Table of attributes associated with each annotation
  %   numAttributes: Number of attributes
  %   attributeNames: Names of attributes
  %
  % Methods:
  %   MSAnnotationSet: Constructor
  %   names: Get annotation names
  %   byName: Get index of an annotation specified by name
  %   find: Find annotations by name matching a regular expression pattern
  %   list: Return table of annotations with index, name and size
  %   items: Get item subsets as logical index mask
  %   segmentItems: Get one or more segments as a logical index mask
  %   segmentLabels: Get label vector for mutually disjoint segments
  %   append: Append annotations
  %   delete: Delete annotations
  %   reduce: Reduce annotation set to specified item subset
  %   assert: Check object consistency, abort if inconsistent
  %
  % MSAnnotationSet uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = immutable)
    numItems  % Total number of items in the dataset
  end
  properties (SetAccess = private)
    annotations = ...  % List of annotations with name and item index set
      struct('name',{},'mask',{});
    attributes = table; % Table of attributes associated with annotations
  end
  properties (Dependent)
    numAnnotations % Number of annotations
    numAttributes  % Number of attributes
    attributeNames % Names of attributes
  end
  
  methods
    function obj = MSAnnotationSet (X)
      % Constructor
      % obj = MSAnnotationSet(numItems): Create empty annotation set for a
      %   mass spec dataset with numItems data items. 
      % obj = MSAnnotationSet(slDumpObject): Initialize annotation set from
      %   slDump object
      
      if isScalarGE1(X)
        % Create empty annotation set with given number of data items
        obj.numItems = floor(X);
      elseif isa(X, 'slDump')
        % Initialize annotation set from slDump object
        obj.numItems = X.getNumberOfSpectra();
        obj.initFromSLObject(X);
      else
        error('Argument must be a positive integer scalar or an slDump object');
      end
      obj.assert();
    end
    
    function N = get.numAnnotations (obj)
      N = length(obj.annotations);
    end
    
    function N = get.numAttributes (obj)
      N = size(obj.attributes, 2);
    end
    
    function N = get.attributeNames (obj)
      N = obj.attributes.Properties.VariableNames;
    end
    
    function N = names (obj, A)
      % Get annotation names
      % N = obj.names(A): Get names of specified annotations as cell string
      %   array. A may be either an index list or a regular expression 
      %   pattern. If A is omitted, all names are returned.
      
      narginchk(1,2);
      if nargin < 2
        A = 1:obj.numAnnotations;
      else
        A = obj.resolveAnnotations(A);
      end
      N = {obj.annotations(A(:)).name}';
    end
    
    function I = byName(obj, name)
      % Get index of an annotation specified by name
      % I = obj.byName(name): Return indices of all annotations with the
      %   exact given name

      if ~ischar(name)
        error('Argument must be a string');
      end
      pattern = strcat('(|/)', regexptranslate('escape', name), '$');
      I = obj.find(pattern);
    end

    function varargout = find(obj, pattern)
      % Find annotations by name matching a regular expression pattern
      % I = obj.find(pattern): Return indices of annotations matching
      %   regexp pattern
      % [I,N] = obj.find(pattern): Return indices and full names

      nargoutchk(0,2);
      if ~ischar(pattern)
        error('Argument must be a string');
      end
      % Find indices of annotation names matching namePattern
      findFn = @(s) ~isempty(regexp(char(s), pattern, 'once'));
      varargout{1} = find(cellfun(findFn, {obj.annotations.name}));
      % Store matching names as optional second output argument
      if nargout >= 2
        varargout{2} = {obj.annotations(varargout{1}).name};
      end
    end

    function T = list(obj, A)
      % Return table of annotations with index, name and size
      % T = obj.list(A): List all specified annotations whith name and size.
      %   A may be either an index list or a regular expression pattern. 
      %   If A is omitted, all annotations are listed.
      
      narginchk(1,2);
      if nargin < 2
        A = 1:obj.numAnnotations;
      else
        A = obj.resolveAnnotations(A);
      end
      names = {obj.annotations(A).name};
      sizes = cellfun(@(M) sum(M>0), {obj.annotations(A).mask});
      T = table(A(:), names', sizes', 'VariableNames', {'Index', 'Name', 'Size'});
    end
    
    function M = items (obj, A)
      % Get item subsets as logical index mask
      % M = obj.items(A): Get items of specified annotations as a logical
      %   m-by-n-matrix, where m is the total number of items (numItems)
      %   and n is the number of annotations. Each column in the result
      %   matrix represents an annotation as a mask vector. 
      %   If A is omitted, item subsets of all annotations are returned.
      
      narginchk(1,2);
      if nargin < 2
        A = 1:obj.numAnnotations;
      else
        A = obj.resolveAnnotations(A);
      end
      M = false(obj.numItems, numel(A));
      for k = 1:numel(A)
        M(:, k) = obj.annotations(A(k)).mask;
      end
    end
    
    function M = segmentItems(obj, varargin)
      % Get one or more segments (sets of annotations) as an index mask
      % M = obj.segmentItems(A): Return a single logical column vector for 
      %   annotations A, specified as an index list or a name pattern.
      %   obj.segmentItems(A) is equivalent to any(obj.items(A),2).
      % M = obj.segmentItems(A1, A2, ...):
      %   Return a logical matrix, where the i-th column represents the 
      %   items in Ai as a logical column vector.
      % M = obj.segmentItems(C) with a cell array vector C is equivalent to
      %   passing the cell array elements as separate arguments.
      
      % Check input arguments
      narginchk(2,inf);
      % Convert input arguments to cell array
      if nargin == 2 && iscell(varargin{1})
        % Single input argument is a cell array, must be a vector
        if ~isvector(varargin{1})
          error('Cell array argument must be 1-by-n or n-by-1 vector')
        else
          C = varargin{1};
        end
      else
        % Multiple arguments
        C = varargin;
      end
      
      % Initialize output matrix
      M = false(obj.numItems, length(C));
      % Iterate over segments
      for k = 1:length(C)
        % Get data item mask for k-th segment
        ind = obj.resolveAnnotations(C{k});
        M(:,k) = any(obj.items(ind),2);
      end
    end

    function labels = segmentLabels(obj, varargin)
      % Get label vector for mutually disjoint segments
      % labels = obj.segmentLabels(A)
      % labels = obj.segmentLabels(A1, A2, ...)
      % labels = obj.segmentLabels(C):
      %   Arguments A, A1, A2, ... or the cell elements of C specifiy
      %   segments (sets of annotations) as with the segmentItems() method.
      %   All forms return an integer column vector of length numItems where
      %   entry values 1, 2, ... correspond to items in segment A1, A2, ...
      %   A value of 0 indicates items that are not in any segment.
      
      M = obj.segmentItems(varargin{:});
      % Check whether segments are disjoint
      if any(sum(M,2) > 1)
        error('Specified annotations are not disjoint')
      end
      % Set output label to k for items in segment k
      labels = zeros(size(M,1),1);
      for k = 1:size(M,2)
        labels(M(:,k)) = k;
      end
    end
    
    function append (obj, varargin)
      % Append annotations
      % obj.append(name, M): Append single annotation
      %   name:  Annotation name (character array)
      %   M: Item subset as a logical vector of length numItems
      % obj.append(names, M): Append multiple annotations
      %   names: Cell array vector of annotation names
      %   M: 2D logical array with numItems rows where the k-th column
      %   corresponds to the k-th name in names
      % obj.append(name(s), M, att): Append annotation(s) with attributes
      %   att, specified as a table object.
      % obj.append(A): Append annotations from MSAnnotationSet A
      
      narginchk(2,4);
      if nargin == 2
        % Single input argument is an MSAnnotationSet object
        if ~isa(varargin{1}, 'MSAnnotationSet')
          error('Argument is not an MSAnnotationSet object');
        else
          A = varargin{1};
          if A.numItems ~= obj.numItems
            error('Number of items do not match');
          else
            n = length(obj.annotations);
            k = length(A.annotations);
            obj.annotations(n+1:n+k,1) = A.annotations;
            if isempty(A.attributes)
              if ~isempty(obj.attributes)
                obj.attributes(n+1:n+k,:) = cell(k, obj.numAttributes);
              end
            else
              obj.attributes(n+1:n+k, A.attributeNames) = A.attributes;
            end
          end
        end
      elseif nargin >= 3
        names = varargin{1};
        M = varargin{2};
        if nargin == 4
          att = varargin{3};
        else
          att = table;
        end
        if ischar(names)
          % Input arguments are a single pair of name and item subset.
          % Convert names to cell string array and mask to column vector.
          names = {names};
          if isvector(M)
            M = M(:);
          end
        end
        if iscellstr(names)
          % Input arguments are a cell array of names and a mask matrix
          if ~isvector(names)
            error('First argument must be a cell string vector');
          elseif ~(islogical(M) && size(M,1) == obj.numItems)
            error(['Second argument must be a logical array with', ...
                   'numItems (= %d) rows'], obj.numItems);
          elseif length(names) ~= size(M,2)
            error('Number of names and columns of mask matrix must match');
          elseif ~(istable(att) && (isempty(att) || size(att,1) == length(names)))
            error(['Third argument must be a table with row count matching ' ...
                   'number of names']);
          else
            if isempty(att)
              att.Variables = cell(length(names),0);
            end
            for k = 1:length(names)
              obj.add(names{k}, M(:,k), att(k,:));
            end
          end
        else
          error('First argument must be a string or a cell array of strings');
        end
      end
    end
    
    function delete (obj, A)
      % Delete annotations
      % obj.delete(A): Delete specified annotations given as indices or as
      %   a regular expression name pattern
      ind = obj.resolveAnnotations(A);
      obj.annotations(ind) = [];
      if ~isempty(obj.attributes)
        obj.attributes(ind,:) = [];
      end
    end
    
    function A = reduce (obj, itemMask)
      % Reduce annotation set to specified item subset
      % A = obj.reduce(itemMask): Reduce annotation set to the subset
      %   specified by the logical or numeric vector itemMask
      
      % Check input argument
      if ~((isnumeric(itemMask) || islogical(itemMask)) && ...
           isvector(itemMask) && length(itemMask) == obj.numItems)
        error(['itemMask must be a numeric or logical vector of length ' ...
               'obj.numItems (=%d)'], obj.numItems);
      end

      % Create a reduced copy of the boolean annotation mask
      annotationMask = obj.items();
      annotationMask = annotationMask(itemMask,:);
      % Remove empty annotations
      nonEmptyAnnotations = find(any(annotationMask,1));
      annotationMask = annotationMask(:, nonEmptyAnnotations);
      annotationNames = obj.names(nonEmptyAnnotations);
      if isempty(obj.attributes)
        att = table;
      else
        att = obj.attributes(nonEmptyAnnotations,:);
      end
      % Create new annotation set object
      A = MSAnnotationSet(sum(itemMask));
      A.append(annotationNames, annotationMask, att);
    end
    
    function assert (obj)
      % Check object consistency, abort if inconsistent
      
      % Check member attributes
      assert(isempty(obj.attributes) || size(obj.attributes,1) == obj.numAnnotations, ...
             ['Assertion failed: Row count of MSAnnotationSet.attributes (%d) ', ...
              'must match number of annotations (%d)'], ...
             size(obj.attributes,1), obj.numAnnotations);
    end
  end
  
  methods (Access = protected)
    function ind = resolveAnnotations (obj, A)
      % Check whether input argument A is an annotation index vector or a
      % string. If it is a string, use it as a name pattern to lookup
      % corresponding annotation indices.
      if ischar(A)
        ind = obj.find(A);
      elseif isIndexVector(A, 1, obj.numAnnotations)
        ind = A;
      else
        error(['Annotation argument must be either a regular expression ', ...
               'name pattern or an integer index vector in range ', ...
               '1..numAnnotations (= %d)'], obj.numAnnotations);
      end
    end
    
    function add (obj, name, M, att)
      % Add single annotation (w/o consistency checks, M is 
      % expected to be a mask vector), att a single row table
      
      n = length(obj.annotations)+1;
      obj.annotations(n,1).name = name;
      obj.annotations(n,1).mask = M(:);
      if ~(isempty(obj.attributes) && isempty(att))
        % Not both attribute tables are empty, append attribute values
        obj.attributes(n, att.Properties.VariableNames) = att;
      end
    end
    
    function initFromSLObject (obj, S)
      % Initialize annotation set from slDump object
      
      for i = 1:length(S.RegionNames)
        M = false(obj.numItems,1);
        M(S.RegionSpots{i}) = true;
        obj.append(S.RegionNames{i}, M);
      end
      obj.attributes = cell2table(S.RegionProperties, ...
                                  'VariableNames', S.RegionPropertyNames);
    end
    
  end
  
end


function result = isScalarGE1 (X)
  result = isnumeric(X) && isscalar(X) && X >= 1;
end

function result = isIndexVector(X, minval, maxval)
  result = isnumeric(X) && isvector(X) && all(X >= minval) && all(X <= maxval);
end
