classdef MSLabelData < MSPredictionData
  % Store label datasets
  %
  % This class stores label assignment data for mass spec datasets. It is
  % derived from the base class MSData and inherits properties for storing
  % data items (label values) and optional related information (annotations, 
  % positions). Label values are assumed to be positive integers where each
  % distinct value is related to an individual label.
  %
  % Properties (in addition to superclass):
  %   labels: Names associated with distinct label values
  %   numLabels: Number of labels
  %   usedLabels: actually occurring labels as a vector of integers 
  %  (labels {k} where obj.data == k is non-empty; found via unique(obj.data))
  %
  % Methods:
  %   MSLabelData: Constructor
  %   items: Get item subsets for each label value
  %   complement: Create complement label assignment
  %   randomSubsets: Generate random, disjoint subsets of specified sizes
  %   kFold: Create label data for k-fold cross validation
  %   leaveOut: Create label data for leave-out cross validation segments
  %   spread: Spread labels to multiple columns
  %   aggregate: Aggregate labels within segments
  %   removeUnused: Remove unused labels
  %
  % Static methods:
  %   getLabelVector: Check and convert argument to column label vector
  %
  % MSLabelData uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.

  properties (Dependent)
    usedLabels % Actually occuring labels (vector of integers), including 0
    usedLabelsExcl0 % Skipping the 0 entry in case it occurs
    numLabelsOccIncl0
    numLabelsOccExcl0
  end
  
  methods
    function obj = MSLabelData (labels, data, sourceData)
      % Constructor
      % obj = MSLabelData(labels, data, sourceData)
      %   labels: Cell string array with names associated with individual 
      %     label values
      %   data: Label data matrix containing label assignments as columns
      %   sourceData (optional): Source object, derived from MSData. If 
      %     specified, positions are copied from sourceData.
      
      narginchk(2,3);
      if nargin < 3
          sourceData = [];
      end
      % Check input arguments
      if ~isempty(sourceData) && ~isa(sourceData, 'MSData')
        error('sourceData must be either empty or an MSData object');
      end
    
      % Call superclass constructor
      obj@MSPredictionData(labels, data);
      % Set labels
      obj.setLabels(labels);
      % If source data is given, copy positions. (Don't copy annotations, 
      % as these are typically not used on derived data objects.)v
      if ~isempty(sourceData) 
        if ~isempty(sourceData.positions)
          obj.setPositions(sourceData.positions);
        end
      end
    end
    
    function usedL = get.usedLabels(obj)
      % Get actually occurring labels
      usedL = unique(obj.data);
    end
    
    function usedLExcl0 = get.usedLabelsExcl0(obj)
      % Get occurring labels but exclude zero in case it occurs
      usedLExcl0 = unique(obj.data);
      usedLExcl0( usedLExcl0 == 0 ) = [];
    end
    
    function n = get.numLabelsOccIncl0(obj)
      n = numel( obj.usedLabels );
    end
    
    function n = get.numLabelsOccExcl0(obj)
      n = numel( obj.usedLabelsExcl0 );
    end
      
      
    function M = items (obj)
      % Get item subsets for each label value as logical index mask
      % M = obj.items(): Get items corresponding to each of the label
      %   values as a logical m-by-n-matrix, where m is the total number of
      %   items (numItems) and n is the number of labels (numLabels). Each
      %   column in the result matrix represents a label segment as a mask
      %   vector.
      
      if obj.dataLength > 1
        error('item() method is not supported for multicolumn labels');
      end
      ind = find(obj.data > 0);
      M = accumarray([ind obj.data(ind)], true, ...
                     [obj.numItems obj.numLabels], @any, false);
    end

    function outLabels = complement (obj, U)
      % Create labels assigned to the complement data item set
      % outLabels = obj.complement(U): Create new MSLabelData object with 
      %   labels assigned from U to the items in obj's complement set.
      %   U may be specified as a column vector or as an MSLabelData object
      %   with one label column. In both cases, number of items must match.
      
      % Check input argument
      if isa(U, 'MSLabelData')
        % Argument is an MSLabelData object
        if ~(U.numItems == obj.numItems && U.dataLength == 1)
          error(['MSLabelData argument must have one label column and the ' ...
                 'same numItems (= %d)'], obj.numItems)
        end
        % Use first data column as basic label set
        U = U.data(:,1);
      elseif isnumeric(U)
        if ~(~U.isempty() && U.isvector() && length(U) == obj.numItems)
          error(['Array argument must be a numeric vector with length ' ...
                 'numItems (= %d)'], obj.numItems)
        end
      else
        error('Argument must be an MSLabelData object or a numeric vector')
      end
      
      % If obj has more than one column, replicate basic labels
      if obj.dataLength > 1
        U = repmat(U, 1, obj.dataLength);
      end
      % Compute complement
      C = zeros(obj.numItems,obj.dataLength);
      mask = obj.data == 0;
      C(mask) = U(mask);
      % Create output label data
      outLabels = MSLabelData(obj.labels, C, obj);
    end
    
    function outLabels = randomSubsets (obj, Nsub, F)
      % Generate random, disjoint subsets of specified sizes
      % outLabels = obj.randomSubsets(Nsub): Create new MSLabelData objects
      %   with labels 1..K defining subsets of sizes as specified in the
      %   vector Nsub (K=length(Nsub)). Subsets are randomly chosen such
      %   that each subset roughly has the same proportion of classes as
      %   represented by the labels in obj. Items without labels are
      %   ignored.
      % outLabels = obj.randomSubsets(Nsub, F): As above, with the vector F
      %   specifying the relative proportions of subset group sizes.

      % Check input arguments
      obj.assert;
      if size(obj.data,2) ~= 1
        error('data must be a single column of label values')
      elseif nargin < 3
        F = [];
      end
      outLabels = MSLabelData(strtrim(cellstr(num2str((1:length(Nsub))'))), ...
                              zeros(obj.numItems,1));
      outLabels.data = MSRandomSubsets(obj.data, Nsub, F);
    end
    
    function outLabels = kFold (obj, K, L)
      % Create label data for k-fold cross validation segments
      % outLabels = obj.kFold (K): Create new MSLabelData object with
      %   labels 1..K defining segments to be used for stratified K-fold
      %   cross validation. Segments are randomly chosen such that each
      %   segment roughly has the same proportion of classes as represented
      %   by the labels in obj.
      % outLabels = obj.kFold (K,L): Align cross validation segments to
      %   partition represented by MSLabelData object L. This is similar to
      %   obj.kFold(K), with the additional restriction that all items with
      %   the same labels in L are assigned to the same cross validation
      %   segment.
      
      % Check input arguments
      obj.assert;
      if size(obj.data,2) ~= 1
        error('data must be a single column of label values')
      elseif ~(isnumeric(K) && isscalar(K) && K >= 1)
        error('k must be a scalar, positive integer')
      elseif nargin < 3
        L = [];
      elseif ~(isa(L, 'MSLabelData') && L.numItems == obj.numItems && ...
           L.dataLength == 1)
        error(['Argument L must be an MSLabelData object with one label ' ...
               'column and the same numItems (= %d)'], obj.numItems);
      elseif any(obj.data > 0 & ~L.data)
        error('Alignment label in L missing for some labeled data items');
      end
      if ~isempty(L)
        L.assert;
      end
      
      % Select non-zero label values
      dataMask = obj.data(:,1) > 0;
      D = obj.data(dataMask,1);

      % Initialize labels vector
      kFoldLabels = zeros(obj.numItems,1);
      if isempty(L)
        % If L was omitted, generate stratified k-fold partition
        kFoldLabels(dataMask,1) = obj.makeKFoldLabels(D,K);
      else
        % Otherwise generate aligned k-fold partition
        kFoldLabels(dataMask,1) = ...
          obj.makeKFoldLabelsAligned(D, K, L.data(dataMask,1), L.numLabels);
      end
      % Create output label data object
      cvLabels = strsplit(strtrim(sprintf('CV-%0*d\n', [ones(1,K)*floor(log10(K))+1; 1:K])));
      outLabels = MSLabelData(cvLabels, kFoldLabels, obj);
    end
    
    function outLabels = leaveOut (obj, L)
      % Create label data for leave-out cross validation segments
      % outLabels = obj.leaveOut(L): Create new MSLabelData object with
      %   segments to be used for leave-out cross validation. Leave-out
      %   subsets are all non-empty intersections of labels in obj with
      %   unique label values in L.
      %   Resulting label values are consecutive numbers 1..(number of
      %   segments). Label names are copied from L.
      
      % Check input arguments
      obj.assert;
      if size(obj.data,2) ~= 1
        error('data must be a single column of label values')
      elseif ~(isa(L, 'MSLabelData') && L.numItems == obj.numItems && ...
           L.dataLength == 1)
        error(['Argument L must be an MSLabelData object with one label ' ...
               'column and the same numItems (= %d)'], obj.numItems);
      elseif any(obj.data > 0 & ~L.data)
        error('Alignment label in L missing for some labeled data items');
      end
      L.assert;
      
      % Limit to data items with non-zero label values
      dataMask = obj.data(:,1) > 0;
      % Find unique label values and index map
      [uniqueValues,~,uniqueMap] = unique(L.data(dataMask,1));
      % Set new labels vector and label names
      leaveOutLabels = zeros(obj.numItems,1);
      leaveOutLabels(dataMask,1) = uniqueMap;
      newLabels = L.labels(uniqueValues);
      outLabels = MSLabelData(newLabels, leaveOutLabels, obj);
    end
    
    function outLabels = spread (obj, L)
      % Spread labels to multiple columns as indicated by labels argument L
      % outLabels = obj.spread(L): Create new MSLabelData object with
      %   as many label columns as the number of labels of the MSLabelData 
      %   input argument L. Label columns in outLabels are generated from
      %   obj's label values, column indices are taken from the values in L.
      
      % Check input arguments
      obj.assert;
      if obj.dataLength ~= 1
        error('Must not have more than one label column')
      elseif ~(isa(L, 'MSLabelData') && L.numItems == obj.numItems && ...
           L.dataLength == 1)
         error(['Argument must be an MSLabelData object with one label ' ...
                'column and the same numItems (= %d)'], obj.numItems);
      end
      L.assert;
      
      % Create output array
      spreadLabels = zeros(obj.numItems, L.numLabels);
      for k = 1:L.numLabels
        % Select items with L label == k
        mask = (L.data == k);
        spreadLabels(mask,k) = obj.data(mask);
      end
      % Create output label data object
      outLabels = MSLabelData(obj.labels, spreadLabels, obj);      
    end
    
    
    function outLabels = aggregate (obj, S, aggregation)
      % Aggregate labels within segments
      % outLabels = obj.aggregate(S): Create new MSLabelData object
      %   representing label values aggregated over all items with the same
      %   label value in segment label object S (MSLabelData).
      %
      %   By default, aggregation is done by selecting the postive label
      %   value that occurs most often among all values with the same S
      %   label value ('max' aggregation). In case more than one label
      %   value occurs the maximum number of times, the lowest value is
      %   selected.
      %
      %   The returned outLabels is an MSLabelData object with label values 
      %   sorted by the corresponding S label value and with 
      %   outLabels.numItems == maximum label value in S.
      %   Items with S label value == 0 are ignored.
      %
      % outLabels = obj.aggregate(S, aggregation): Use specified
      %   aggregation method to aggregate labels over segments. aggregation
      %   may be one of the folowing:
      %     'max' - Select positive label value occuring most often within
      %             each segment. If not unique, select the lowest such
      %             label. Items with label value == 0 are ignored.
      %     'unique' - Select 0 if all label values within segment equal 0,
      %                select NaN if >= 2 different values > 0 occur within segment
      %                  (this raises a warning),
      %                otherwise select the unique positive value occuring in segment.
      %     fn - A function handle fn that is called for each segment with
      %          two arguments. The first argument is a vector containing
      %          the number of occurences of all possible positive label 
      %          values, the second argument is the total number of items
      %          in the segment (including items with label value == 0).
      %          fn must return the selected label value as a scalar.
      
      % Check arguments
      obj.assert;
      narginchk(1,3);
      if nargin < 2 || isempty( S )
        S = obj.trivialAggregationLabelData;
      end
      if ~isa(S, 'MSLabelData')
        error('Argument must be an MSLabelData object');
      end
      if nargin < 3
        aggregation = [];
      end
      if isempty(aggregation)
        % Use max aggregation by default
        aggregation = 'max';
      end
      S.assert;
      S_numLabels = S.numLabels;
      S = obj.checkLabelArgument(S);
      
      % Initialize output label object
      outLabels = MSLabelData(obj.labels, zeros(S_numLabels, obj.dataLength));
      % Iterate over label columns
      for col = 1:obj.dataLength
        % Compute contigency table, first column counts label values == 0
        itemMask = S(:, col) > 0;
        C = accumarray([S(itemMask, col), obj.data(itemMask, col)+1], 1, ...
                       [S_numLabels, obj.numLabels+1]);

        % Find segments with non-zero label values
        segmentMask = any(C(:,2:end),2);
        % Apply specified aggregation method
        if ischar(aggregation)
          if strcmp(aggregation, 'max')
            % 'max': Select label (= column index) with highest count
            [~, outLabels.data(segmentMask, col)] = max(C(segmentMask, 2:end), [], 2);

          elseif strcmp(aggregation, 'unique')
            % 'unique': Select unique label, assign NaN if not unique
            [maxCount,maxLabel] = max(C(segmentMask, 2:end), [], 2);
            nonUnique = (sum(C(segmentMask, 2:end), 2) - maxCount) > 0;
            if any(nonUnique)
              maxLabel(nonUnique) = NaN;
              warning(['Non-unique labels found in %d of %d segments, '...
                       'aggregation set to NaN'], sum(nonUnique), length(nonUnique));
            end
            outLabels.data(segmentMask, col) = maxLabel;
          else
            error('Unknown aggregation method');
          end

        elseif isa(aggregation, 'function_handle')
          % Aggregation argument is a function handle, call on each segment
          for k = find(segmentMask)'
            % Call with label counts and total item count
            outLabels.data(k, col) = aggregation(C(k,2:end), sum(C(k,:)));
          end
        elseif isa(aggregation, 'MSAggregationSolver')
         % aggregation is an aggregation solver object. Call resolveAggregation
         % on each segment
          for k = find(segmentMask)'
             outLabels.data(k,col) = ...
                   aggregation.resolveAggregation(C(k,2:end), sum(C(k,:)));
          end
        else
            error('Unsupported aggregation argument type')
        end
      end
    end
    
    function C = confusionMatrix(obj, groundTruth, data)        
        % Specify explicit label value order, including zero values, to
        % make sure that final row/column indices represent label values.
        C = accumarray([groundTruth(:) data(:)]+1, 1, [1 1]*(obj.numLabels+1));
        % Ignore zero items by removing first column and row
        C=C(2:end,2:end);
    end
    
    function varargout = removeUnused (obj)
      % Remove unused labels
      % obj.removeUnused(): Modify obj by removing all label values that
      %   are not used. The relative ordering of the remaining label values
      %   preserved, but values may change.
      % L = obj.removeUnused(): Return new label object with unused label
      %   values removed, leave obj unchanged.
      
      % Assert consistency
      obj.assert;
      % Check number of output arguments
      nargoutchk(0,1);
      % List of used, non-zero label values
      indUsed = obj.usedLabels(:)';
      indUsed = indUsed(indUsed > 0);
      % Label value map converting old to new label values
      % First entry corresponds to label value == 0 (no label)
      labelMap = nan([1 obj.numLabels+1]);
      labelMap(1) = 0;
      labelMap(indUsed+1) = 1:length(indUsed);
      % Apply label map
      D = reshape(labelMap(obj.data+1), size(obj.data));
      % Remove unused labels
      L = obj.labels(indUsed);

      if nargout >= 1
        % Create new output object
        varargout{1} = MSLabelData(L, D, obj);
      else
        obj.data = D;
        obj.labels = L;
      end
    end
    
    function assert (obj)
      % Check object consistency, abort if inconsistent
      
      % Call superclass assert
      assert@MSPredictionData(obj);
      
      % Data values must be integers in 0..numLabels
      usedL = obj.usedLabels;
      assert(all(usedL >= 0 & usedL <= obj.numLabels & mod(usedL,1) == 0), ...
             ['Assertion failed: Label data values must be integers in ' ...
              '0..numLabels (=%d)'], obj.numLabels);
    end
  end
  
  methods (Access = protected)  
    
    function indexVector=vectorFromIndexes(~,index)
        indexVector=index;
    end
    
    function value=getNumPredictions(obj)
        value=obj.dataLength;
    end
    
    function C = getColorMap (obj)
      % Return the color map to be used by the show method
      
      % Use the linspecer colormap with 
      % unused labels (value == 0) set to light gray
      C = [0.77*ones(1,3); linspecer(obj.numLabels)];
    end
    
    function showImage (obj, X)
      % Show vector X as image
      showImage@MSData(obj, X); % call superclass method
      % Make sure that label values map to color map entries, label value 0
      % maps to background color
      caxis([0 obj.numLabels]);
      % Create a legend if there are less than 25 distinct label values
      if obj.numLabels < 25
        hold on;
        N = obj.numLabels;
        % find any point in the foreground
        [rndPtx,rndPty] = ind2sub(obj.positions.gridSize, obj.positions.reverseIndex(1));
        % create N dummy lines from the rndPoint to itself
        L = line(repmat(rndPtx,2,N),repmat(rndPty,2,N),'LineWidth',2);
        cmap = colormap(gca);
        set(L,{'color'},mat2cell(cmap(2:(N+1),:),ones(1,N),3));
        % set the colors according to cmap (cmap in MSLabelData is used, first entry reserved
        % for background)
        legend(L,obj.labels,'location','southeast','interpreter','none'); % add legend for the dummy lines
        hold off
      end
    end
    
    function outLabels = select_impl(obj, N)
      % Select a subset of label columns from a label matrix
      % outLabels = obj.select(N): Return a new MSLabelData object
      %   representing the N-th label column(s) of obj     
      outLabels = MSLabelData(obj.labels, obj.data(:,N(:)), obj);
    end
    
    function D = newFromTemplate (obj, data)
      % Create a new MSLabelData object with specified data
      D = MSLabelData(obj.labels, data);
    end
    
    function updateLabels_impl (obj, fwdLabelMap, ~)
      % Apply mapping to data values, unless mapping vector equals (1:n)
      if any(fwdLabelMap ~= (1:length(fwdLabelMap))')
        fwdLabelMap = [0; fwdLabelMap];
        obj.data = fwdLabelMap(obj.data+1);
      end
    end
  end
  
  methods (Static)
    function [L,N,usedLabels] = getLabelVector (labels, msData)
      % Check and convert argument to column label vector
      % L = getLabelVector(labels): Argument labels may be an MSLabelData
      %   object with a single label column or a numeric vector.
      % L = getLabelVector(labels, msData): If the optional msData
      %   argument is given, length of labels is checked to match
      %   msData.numItems.
      % [L,N] = getLabelVector(...): Additionally return number of labels
      % [L,N,usedLabels] = getLabelVector(...): Additionally return number of
      %     labels and the actually used labels 
      % (labels {k} where obj.data == k is non-empty; found via unique(obj.data))
      
      narginchk(1,2)
      nargoutchk(0,3)
      if isa(labels, 'MSLabelData')
        labels.assert;
        if labels.dataLength ~= 1
          error('MSLabelData argument must have one label column')
        end
        L = labels.data(:,1);
      elseif isnumeric(labels) && ~isempty(labels) && isvector(labels)
        L = labels(:);
      else
        error(['labels argument must be an MSLabelData object or a', ... 
            'non-empty numeric vector']);
      end
      if nargin >= 2
        if ~isa(msData, 'MSData')
          error('msData argument must be an MSData object');
        elseif msData.numItems ~= size(L,1)
          error('Length of label vector must match number of items in msData argument')
        end
      end
      % If requested, return number of labels in second output argument
      if nargout >= 2
        if isa(labels, 'MSLabelData')
          N = labels.numLabels;
        else
          N = max(floor(max(L)),0);
        end
      end
      % If requested, return usedLabels in third output argument
      if nargout >= 3
        if isa(labels, 'MSLabelData')
          usedLabels = labels.usedLabels;
        else
          usedLabels = unique(L);
        end
      end
    end
  end
  
  
  methods (Static, Access = protected)
    function labels = makeKFoldLabels (D, K)
      % Make label vector for stratified k-fold cross validation.
      % D is expected to be a column vector of positive integers
      % representing the stratification groups. K is a positive integer.
      
      % Create k-fold partition
      CV = cvpartition(D, 'KFold', K);
      % For each of the K repetitions, assign repetition number as label 
      % value to the items in the corresponding test set.
      labels = zeros(size(D));
      for i = 1:CV.NumTestSets
        labels(CV.test(i),1) = i;
      end
    end
    
    function labels = makeKFoldLabelsAligned (D, K, L, Lmax)
      % Make label vector for aligned, stratified k-fold cross validation.
      % D is expected to be a column vector of positive integers
      % representing the stratification groups. K is a positive integer.
      % L is a column vector of positive integers in the range 1..Lmax 
      % representing the segments to align with.
      
      % Compute a contingency table where
      %   T(i,j) = Number of items with L = values{i,1} and D = values{j,2}
      [T,~,~,values] = crosstab(L,D);
      % For each label value in L, find most frequent label in D
      [~,Dmajor] = max(T,[],2);
      % Generate stratified k-fold cross validation partition on L-segments
      segmentLabels = MSLabelData.makeKFoldLabels(Dmajor,K);
      % Some L-segments may not be included in the CV partition, as they
      % may have no intersection with the labels in D. We need to map
      % between CV indices and L labels.
      LToCVMap = zeros(Lmax,1);
      for i = 1:size(values,1)
        if ~isempty(values{i,1})
          LToCVMap(str2double(values{i,1})) = i;
        end
      end
      % Map segment labels to individual item labels
      labels = segmentLabels(LToCVMap(L));
    end
  end
end

