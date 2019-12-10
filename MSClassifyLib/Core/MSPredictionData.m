classdef MSPredictionData < MSData
    % This class (abstract) represents a classification prediction. So far
    % it can be either a deterministic classification (subclass
    % MSLabelData) or a fuzzy classification (subclass MSScoreData)
    %
    % Properties
    %   labels: (MSLabelData representing the classification classes or
    %           labels)
    %   numLabels: Number of classes or labels
    %   numPredictions: Number of different predictions stored
    %
    % Methods:
    %   MSPredictionData: Constructor
    %   setlabels: Sets labels object
    %   aggregateInverse: Expands prediction over segments (revert aggregation)
    %   chekLabelArgument: Validate label for aggregation
    %   getLabelMatrix: Convert label into matrix for aggregation
    %   validateIndex: Produces an error if the input index is not a valid
    %        prediction index vector (used for example in select)
    %   select: Creates an MSPredictionData from the data indicated by the
    %       input index
    %   select_impl (abstract): Creates the MSPredictionData corresponding
    %       to the output of 'select'
    %   aggregate (abstract): Aggregates predictions within segments
    %   getNumPredictions (abstract): Gets number of predictions
    %   join (static): Combine cell array of predictions into a single object

    
    properties (SetAccess=protected)
        labels %(MSLabelData representing the classification classes or
    %           labels)
    end
    properties (Dependent)
        numLabels % Number of classes or labels
        numPredictions % Number of different predictions stored
        anyPredictionsMask % Mask of rows for which predictions (non-zero) 
        % exist
    end
    
    methods
      function obj = MSPredictionData (labels, data, sourceData)
      % Constructor
      % obj = MSPredictionData(labels, data, sourceData)
      %   labels: Cell string array with names associated with individual 
      %     label values
      %   data: Label data matrix containing label assignments as columns
      %   sourceData (optional): Source object, derived from MSData. If 
      %     specified, annotations and positions are copied.
      
      narginchk(2,3);
      if nargin < 3
          sourceData = [];
      end
      % Check input arguments
      if ~isempty(sourceData) && ~isa(sourceData, 'MSData')
        error('sourceData must be either empty or an MSData object');
      end      

      data = single(data); % Reducing memory issues of large 
      % MSClassificationResult objects by 50% (ScoreData).
      % Alternatively this line could be added in the MSScoreData
      % constructor, and MSLabelData can be set to something like uint16.
      
      % Call superclass constructor
      obj@MSData(data);
      % Set labels
      obj.setLabels(labels);
      % If source data is given, copy annotations and positions
      if ~isempty(sourceData) 
        if ~isempty(sourceData.positions)
          obj.setPositions(sourceData.positions);
        end
      end
      end
      
      function N = get.numLabels (obj)
          % Get number of labels
          N = numel(obj.labels);
      end
      
      function value=get.numPredictions(obj)
          % gets number of predictions (depends on the abstract getNumPredictions)
          value=getNumPredictions(obj);
      end
      
      function mask = get.anyPredictionsMask(obj)
        mask = any( logical(obj.data), 2 );
      end
      
      function setLabels (obj, labels)
          % Set names associated with distinct label values

          if ~(iscellstr(labels) && isvector(labels) && ~isempty(labels))
              error('labels must be a non-empty cell string array')
          else
              % Set labels
              obj.labels = labels;
          end
      end  
      
      function updateLabels (obj, newLabels, allowRemove)
          % Update labels and adjust label values accordingly
          % obj.updateLabels(newLabels): Set labels to newLabels, which must be
          %   a superset of previous labels, possibly in modified order. Label
          %   data values are modified accordingly in case a label string is
          %   moved to a different position.
          % obj.updateLabels(newLabels, allowRemove): Same as above, but
          %   specify whether previously existing labels may be omitted in 
          %   newLabels, causing the corresponding data values to be set to 0.
          %   allowRemove may be either of 'allowRemove', true, or false, the
          %   default is false.

          % Check input arguments
          narginchk(2,3);
          if ~(iscellstr(newLabels) && isvector(newLabels) && ~isempty(newLabels))
              error('newLabels must be a non-empty cell string array')
          end
          if nargin < 3
             allowRemove = [];
          elseif ischar(allowRemove)
              if ~(isvector(allowRemove) && strcmpi(allowRemove, 'allowRemove'))
                  error('Invalid allowRemove argument');
              end
              allowRemove = true;
          elseif ~(isempty(allowRemove) || (islogical(allowRemove) && isscalar(allowRemove)))
              error('Invalid allowRemove argument');
          end
          if isempty(allowRemove)
             allowRemove = false;
          end

          if ~allowRemove && ~all(ismember(obj.labels, newLabels))
              error('New labels must contain all current labels');
          end

          % Create matching matrix
          labelMatch = false(obj.numLabels, length(newLabels));
          for k = 1:length(newLabels)
              labelMatch(:,k) = strcmp(obj.labels, newLabels{k});
          end
          % Verify that mapping of current to new labels is unique
          if any(sum(labelMatch,2) > 1)
              error('New labels are not unique');
          end

          % Generate vector mappings between current and new label values
          [iCur, iNew] = find(labelMatch);
          fwdLabelMap = zeros(obj.numLabels,1);
          fwdLabelMap(iCur) = iNew;
          bwdLabelMap = zeros(length(newLabels),1);
          bwdLabelMap(iNew) = iCur;
          
          % Apply mapping to data
          obj.updateLabels_impl(fwdLabelMap, bwdLabelMap);
          % Set new labels
          obj.setLabels(newLabels);
      end
    
      function out = aggregateInverse (obj, S)
      % Revert the aggregate method above in the sense that an MSLabelData/
      % MSScoreData object with the original data length is produced such that all
      % items belonging to one segment all copy the label assignement from
      % the associated aggregation label object:
      %
      % [ 0 4 1 4 0 2 3 3 3 0 ]     (samples = [ 0 1 1 1 0 2 2 2 2 0 ])
      %                 |  
      %                 | (aggregate 'max')
      %                v
      %             [ 4 3 ]                 (samples = [ 0 1 1 1 0 2 2 2 2 0 ])
      %                 |  
      %                 | (aggregateInverse)
      %                v
      % [ 0 4 4 4 0 3 3 3 3 0 ]
      %
      % - if used in the standard way S has to be the same segment label object
      %     previously used for aggregation
      % - by construction aggregateInverse does not depend on the choice of
      %     the used aggregation method ('max', 'unique' or fnhandle), i.e.
      %     there is no method choice aggregationInverse
      % - the case obj.dataLength > 1 is allowed, aggregateInverse will
      %     operate on the columns independently in this case
      %
      % if there is a pixel-wise prediction it is possible to directly
      % take a look at the aggregated prediction as follows:
      %     predictionAgg = prediction.aggregate(S);
      %     predictionAggI = predictionAgg.aggregateInverse(S);
      %     predictionAggI.show;
      %

      obj.assert;
      if ~(isa(S, 'MSLabelData') && S.numLabels == obj.numItems ...
              && S.dataLength == 1)
         error(['Argument S must be a 1-column MSLabelData object with ' ...
                'number of labels equal to numItems (= %d)'], obj.numItems);
      end
      S.assert;

      itemMask = S.data(:,1) > 0;
      L = zeros(S.numItems, obj.dataLength);
      L(itemMask,:)=obj.data(S.data(itemMask),:);
      if isa(obj,'MSLabelData')
        out = MSLabelData(obj.labels, L, S);
      elseif isa(obj,'MSScoreData')
        out = MSScoreData(obj.labels, L, S);
      end
      
      end
      
      function predictionData=select(obj,index)
          obj.validateIndex(index);
          predictionData=obj.select_impl(index);
      end
      
      function data=getPredictionData(obj,index,checkValidity)
            % gets the prediction data for specified index
            % data=obj.getPredictionData(index, checkValidity)
            % INPUT
            %   index: vector with indexes of selected predictions
            %   checkValidity: logical indicating whether the validity of
            %   the index vector should be checked or not (optional, by
            %   default it is true)
            % OUTPUT
            %   data: matrix containing prediction
            
            narginchk(2,3)
            if nargin<3 
                checkValidity=true;
            elseif isempty(checkValidity)||~(isscalar(checkValidity)&& islogical(checkValidity))
                error('Input checkValidity must be a logical scalar')
            end
            if checkValidity
                obj.validateIndex(index);
            end
            % build composed index vector
            indexVector=obj.vectorFromIndexes(index);            
            %return data in selected columns
            data=obj.data(:,indexVector);
      end
      
      function varargout = validate (obj, R)
      % Compute agreement with reference labels
      % P = obj.validate(R): Compute accuracy P of agreement to reference 
      %   labels R as the ratio 
      %     number of equal non-zero labels / number of non-zero labels
      %   R may be one of the following:
      %     - MSLabelData object with one label column
      %     - MSLabelData object with same number of label columns as obj
      %     - Numeric vector with numItems entries
      %     - Numeric m-by-n matrix with m == numItems and n == number of 
      %       label columns
      %
      %   If obj has multiple columns accuracy is returned as a vector with
      %   a separate value for each column. If R has one label column, R is 
      %   considered to be the reference data for each label column of obj. 
      %   If R has multiple columns, each column is considered to be the
      %   reference for the corresponding column of obj.
      %
      % [P, C] = obj.validate(R): In addition to the performance returned 
      %    in P, return the confusion matrix in C, counting the number of
      %    items with a given combination of observed and reference label
      %    value. Rows in C correspond to reference labels (R values), 
      %    columns correspond to observed label values (obj values). If obj
      %    has a single label column, C is a k-by-k square matrix. If obj 
      %    has n label columns, C is a 3D-array of size [k,k,n].
      
      % Check arguments
      obj.assert;
      nargoutchk(0,2)
      R = obj.checkLabelArgument(R);

      C = zeros(obj.numLabels, obj.numLabels, obj.numPredictions);
      accuracy=zeros(obj.numPredictions,1);
      % Confusion matrix is computed independently for each label column
      for k = 1:obj.numPredictions
        C(:,:,k) = obj.confusionMatrix(R(:,k), obj.getPredictionData(k));
        % accuracy is computed
        accuracy(k)=sum(diag(C(:,:,k)))/sum(sum(C(:,:,k)));
      end
      varargout{1}=accuracy;
      if nargout==2
          varargout{2} = C;
      end
      end
      
      function labelData = trivialAggregationLabelData(obj)
        labelData = MSLabelData( {'onlyLabel'}, ones( obj.numItems,1 ), [] );
      end
      
    end
    
    methods(Access = protected)
        function [L,N] = checkLabelArgument (obj, labels)
          % Convert argument to label matrix matching obj's number of labels
          % L = obj.checkLabelArgument(labels): Argument labels may be one of
          %     the following:
          %     - MSLabelData object with one label column
          %     - MSLabelData object with same number of label columns as obj
          %     - Numeric vector with numItems entries
          %     - Numeric m-by-n matrix with m == numItems and n == number of 
          %       label columns
          %     Return value L is an m-by-n matrix with m == numItems and
          %     n == number of label columns. If labels n > 1 and labels has
          %     only one column, its values are replicated.
          % [L,N] = checkLabelArgument(labels): Additionally return number of labels

          nargoutchk(0,2)
          [L,NN] = MSPredictionData.getLabelMatrix(labels, obj);
          if ~((size(L,2) == 1) || (size(L,2) == obj.numPredictions))
            error('Size of labels argument must match number of label columns (= %d)', ...
                  obj.dataLength)
          end
          % If L is a single column and obj has multiple label columns,
          % replicate to match matrix sizes
          if (obj.dataLength > 1) && (size(L,2) == 1)
            L = repmat(L, 1, obj.numPredictions);
          end
          if nargout >= 2
            N = NN;
          end
        end
        
        function validateIndex(obj, index)
            if isempty(index) || ~(isnumeric(index)&& isvector(index)&&...
                    all(index >= 1 & index <= obj.numPredictions))
                error('Argument must be an integer array between 1 and numPredictions (= %d)', ...
                      obj.numPredictions);
            end 
        end
    end
    
    methods(Static)
        function [L,N,usedLabels] = getLabelMatrix (labels, msData)
          % Check and convert argument to label matrix
          % L = getLabelMatrix(labels): Argument labels may be an MSLabelData
          %   object with any number of label columns or a numeric matrix.
          % L = getLabelMatrix(labels, msData): If the optional msData
          %   argument is given, length of labels is checked to match
          %   msData.numItems.
          % [L,N] = getLabelMatrix(...): Additionally return number of labels
          % [L,N,usedLabels] = getLabelMatrix(...): Additionally return number of
          %     labels and the actually used labels
          % (labels {k} where obj.data == k is non-empty; found via unique(obj.data))

          narginchk(1,2)     
          if isa(labels, 'MSLabelData')
            labels.assert;
            L = labels.data;
          elseif isnumeric(labels) && ~isempty(labels) && ismatrix(labels)
            L = labels;
          else
            error(['labels argument must be an MSLabelData object or a ',...
                'non-empty numeric matrix']);
          end
          if nargin >= 2 && msData.numItems ~= size(L,1)
            error('Length of label vector must match number of items in msData argument')
          end
          % If requested, return number of labels in second output argument
          if nargout >= 2
            if isa(labels, 'MSLabelData')
              N = labels.numLabels;
            else
              N = max(floor(max(L(:))),0);
            end
          end
          % If requested, return usedLabels in third output argument
          if nargout >= 3
            if isa(labels, 'MSLabelData')
              usedLabels = labels.usedLabels;
            else
              usedLabels = unique(L(:));
            end
          end
        end
        
        function Q = join (P)
          % Combine cell array of predictions into a single object
          if isa(P, 'MSPredictionData')
            Q = P;
            return;
          end
          if ~(iscell(P) && ~isempty(P) && all(cellfun(@(x) isa(x, 'MSPredictionData'), P)))
            error('Argument must be an array of prediction data objects');
          end
          Q = [];
          for k = 1:numel(P)
            Pk = P{k};
            if isempty(Q)
              Q = P{k}.copy();
            else
              if ~(strcmp(class(Q), class(Pk)) && Q.numItems == Pk.numItems && ...
                   Q.numLabels == Pk.numLabels && Q.dataLength == Pk.dataLength && ...
                   Q.numPredictions == Pk.numPredictions)
                error('Prediction data objects must have same type, size and value range');
              end
              mask = Pk.data > 0;
              if any(Q.data(mask))
                error('Prediction data must be disjoint');
              end
              Q.data(mask) = Pk.data(mask);
            end
          end
        end
    end
    
    methods(Abstract)
        outLabels = aggregate (obj, S, aggregation);
        C = confusionMatrix(obj, groundTruth, data);
    end
    
    methods(Abstract,Access=protected)
        value = getNumPredictions(obj);
        predictionData=select_impl(obj,index);
        indexVector=vectorFromIndexes(obj,index);
        updateLabels_impl(obj, fwdLabelMap, bwdLabelMap);
    end
end