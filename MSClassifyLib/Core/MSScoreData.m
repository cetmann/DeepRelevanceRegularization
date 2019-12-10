classdef MSScoreData < MSPredictionData
    % Class representing a fuzzy classification, where each row of the data
    % object contains the discrete probability distribution of the
    % corresponding spot
    %
    % Properties
    %   numClasses: Number of classes (classification labels)
    %
    % Methods
    %   MSScoreData: Constructor
    %   getPredictionData: gets a section of the data corresponding to the columns
    %       related with a score prediction indicated by an index
    %   setPredictionData: gets a section of the data corresponding to the columns
    %       related with a score prediction indicated by an index
    %   splitScores: returns a cell array of MSScoreData objects with
    %       individual predictions
    %   score2Label: converts the score data into label data
    %   aggregate: aggregate scores within segments
    %   confusionMatrix: compute confusion matrix
    %   ROCCurve: compute ROC curve
    %   ROC_AUC: Compute ROC AUC from scores and reference label L
    %   getNumPredictions: gets number of different predictions stored in
    %       the object
    %   select_impl: implementation of select method (in base class)
    
    properties(Dependent)
        numClasses; %: Number of classes (classification labels)
    end
    methods
        function obj=MSScoreData(labels, data, sourceData)
            % Constructor
            % obj=MSScoreData(labels, data, sourceData)
            % INPUT
            %   labels: cellstr with the names of class labels
            %   data: numerical matrix representing scores. The number of
            %   columns must be a multiple of the number of classes (more
            %   than one prediction can be stored in an object, for
            %   example, useful for storing predictions with different
            %   number of features)
            %   sourceData: optional msdata to copy the positions in order
            %   to be able to show the prediction
            % OUTPUT
            %   obj: MSScoreData
            narginchk(2,3);
            if nargin < 3
                sourceData = [];
            end
            obj@MSPredictionData(labels, data, sourceData);
            obj.assert;
        end

        function setPredictionData(obj,index,data)
            % sets input score data for specified index
            
            %validate index
            obj.validateIndex(index);
            %validate data
            
            %create vector of indexes
            indexVector=obj.vectorFromIndexes(index);
            %assign data
            obj.data(:,indexVector)=data;
        end
        
        function cellScores=splitScores(obj)
            % split object into <numPredictions> objects containing
            % single-score information (for example, a prediction for a
            % corresponding number of features used for training
            cellScores=cell(1,obj.numPredictions);
            for i=1:obj.numPredictions
                cellScores{i}=MSScoreData(obj.labels,obj.getPredictionData(i),obj);
            end
        end 
        
        function labelData=score2Label(obj, fun)
            % Converts the score data into label data
            % labelData=obj.score2Label % uses 'max' probability criterium
            % by default
            % labelData=obj.score2Label(fun) % uses a specified conversion
            % criterium
            % INPUT
            %   fun: Conversion criterium. It can be:
            %       -string:
            %           'max': for each data row the class with the highest
            %               probability is assigned
            %           'prob': for each row a class is asigned according
            %               to the probability distribution represented by 
            %               the row
            %       -function handle:
            %           a function that receives a score matrix and returns
            %           a vector with the class asignation
            
            %check input arguments
            narginchk(1,2)
            if nargin <2
                % by default max criterium is applied
                fun=@obj.score2LabelMax;
            elseif ischar(fun)
                if strcmp(fun,'max')
                    fun=@obj.score2LabelMax;
                elseif strcmp(fun,'prob')
                    fun=@obj.score2LabelProb;
                else
                    error('Unknown conversion type')
                end
            elseif ~isa(fun,'function_handle')
                error('the input <fun> must be a valid function handle')
            end
            % init output labels
            labelData=MSLabelData(obj.labels,zeros(obj.numItems,obj.numPredictions),obj);
            % compute labels for each prediction in obj
            for i=1:obj.numPredictions
                try 
                    labelVector=fun(obj.getPredictionData(i));
                    labelData.data(:,i)=labelVector;
                catch exception
                    msg=exception.message;
                    error(['Wrong function evaluation. Make sure that '...
                        'the input function receives a score matrix '...
                        'and returns a vector where each entry is a '...
                        'number between 0 and the number of classes ('...
                        'numClasses). More about the error: ' msg])
                end                
            end
        end
        
        function value=get.numClasses(obj)
            % gets number of classes
            value=length(obj.labels);
        end
        
        function outScores = aggregate (obj, S, aggregation)
          % Aggregate scores within segments
          % outScores = obj.aggregate(S): Create new MSScoreData object
          %   representing score values aggregated over all items with the same
          %   label value in segment label object S (MSLabelData).
          %
          %   By default, aggregation is done by computing the average
          %   score over all values with the same S label value ('mean' 
          %   aggregation). 
          %
          %   The returned outScores is an MSScoreData object with score values 
          %   sorted by the corresponding S label value and with 
          %   outLabels.numItems == maximum label value in S.
          %   Items with S label value == 0 are ignored.
          %
          % outScores = obj.aggregate(S, aggregation): Use specified
          %   aggregation method to aggregate scores over segments. aggregation
          %   may be one of the folowing:
          %     'mean'   - Compute the average score over each segment
          %     'median' - Compute the median score over each segment
          %
          %     (if it would be useful to have other obvious ones as 'built-in' 
          %     implementedthen other strings can be included in the future.
          %
          %     fn - A function handle fn that is called for each segment with
          %          four arguments:
          %         -data- data score for a single prediction
          %         -S- Label data (single column) representing aggregation
          %             labels
          %         -mask-Logical vector indicating which predictions are
          %             to be taken into account
          %         -numLabels- number of classes represented in S

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
            % Use mean aggregation by default
            aggregation = 'mean';
          end
          S.assert;
          s_numLabels = S.numLabels;
          S = obj.checkLabelArgument(S);
          
          if ischar(aggregation)
              if ~any(strcmp(aggregation, {'mean','median'}))
                  error('Unknown aggregation type')
              end
          elseif ~isa(aggregation, 'function_handle')
              error('Unsupported aggregation argument type')
          end
          %init output score data
          newData=zeros(s_numLabels,obj.dataLength);  
          % Iterate over scores
          for col = 1:obj.numPredictions
            % mask of used labels in aggregation
            itemMask = S(:, col) > 0;  
            % subdata corresponding to current prediction
            data=obj.getPredictionData(col);
            % mask of valid scores (rules out zero rows) 
            maskScore=sum(data,2)~=0;
            % combined mask
            mask=itemMask & maskScore; 
            colIndex=(col-1)*obj.numClasses+1:col*obj.numClasses;
            if ischar(aggregation) % built-in aggregations
                if strcmp(aggregation, 'mean')
                    builtAggreg=@obj.aggregateMean;
                else
                    builtAggreg=@obj.aggregateMedian;
                end
                newData(:,colIndex) = builtAggreg(data,S(:,col),mask,s_numLabels);
            else %user-defined aggregation
                for i=1:s_numLabels
                    itemMask=S(:,col)==i;
                    mask=itemMask & maskScore;
                    if any(mask)
                        newData(i,colIndex)=aggregation(data(mask,colIndex));
                    end
                end
            end
          end   
          outScores=MSScoreData(obj.labels,newData);
        end
        
        function C = confusionMatrix(obj, groundTruth, data)
            % Compute confusion matrix
            
            nclasses=obj.numLabels;
            C = zeros(nclasses);
            for i=1:nclasses
                mask=groundTruth==i;
                C(i,:)=sum(data(mask,:),1);
                % If the dimension 1 is not added but left out since it's
                % the default summing dimension this is different in case
                % mask is a logical vector with 1 true element: then the
                % sum is determined along the 2nd dimension and C(i,:) is
                % filled with constant values! Especially for aggregated
                % data this case is not unlikely.
            end
        end
        
        function R = ROCCurve(obj, L)
            % Compute ROC curve object from scores and reference label L
            % R = obj.ROCCurve(L): Compute MSROCCurve object R for all
            %   predictions in obj, considering the first class as
            %   'positive'. Reference ground truth L may be a logical or
            %   numeric vector or a single column MSLabelData object with the
            %   same number of labels as obj.
            %   Note: Only supported for obj.numLabels == 2
            %

            if obj.numLabels ~= 2
                error('Computing ROC curves for > 2 classes is not supported');
            end
            if isa(L, 'MSLabelData')
                if ~(L.numPredictions == 1 && L.numItems == obj.numItems && ...
                     L.numLabels == obj.numLabels)
                    error(['Label data object L must be single column and must ' ...
                            'match with respect to number of classes and items']);
                end
                L = L.data;
            elseif isnumeric(L)
                if ~(isvector(L) && length(L) == obj.numItems && min(L) >= 0 && ...
                     max(L) <= obj.numLabels && all(mod(L,1) == 0))
                    error(['Numeric label data must be a non-negative integer vector '
                           'and must match with respect to number of classes and items']);
                end
                L = L(:);
            elseif islogical(L)
                if ~(isvector(L) && length(L) == obj.numItems)
                    error('Logical label data must be a vector matching number of items');
                end
                % Convert logical (true,false) to numeric (1,2)
                L = ~L(:)+1;
            end

            % Use only items with reference label and prediction data
            M = any(obj.data > 0, 2) & (L > 0);
            % Compute ROC curves from score data for class 1
            R = MSROCCurve(obj.data(M, 1:2:end), L(M) == 1);
        end
        
        function [Y,A] = ROC_AUC(obj, L)
            % Compute ROC AUC from scores and reference label L
            % Y = obj.ROC_AUC(L): Compute area under curve (AUC) values for
            %   all predictions in obj. Reference ground truth L may be a
            %   numeric vector or a single column MSLabelData object with
            %   the same number of labels as obj.
            % [Y,A] = obj.ROC_AUC(L): Return matrix of pair-wise AUC values
            %   in A (p-by-p-by-n matrix, p = number of classes, n = number
            %   of predictions)
            %
            %   Note: AUC is computed based on a multi-class generalization
            %   of the classical AUC, as described in:
            %   Hand et.al., A Simple Generalisation of the Area Under the
            %   ROC Curve for Multiple Class Classification Problems, 2001.
            %   For two classes, this is equivalent to the classical AUC
            %   definition.
            %

            if isa(L, 'MSLabelData')
                if ~(L.numPredictions == 1 && L.numItems == obj.numItems && ...
                     L.numLabels == obj.numLabels)
                    error(['Label data object L must be single column and must ' ...
                            'match with respect to number of classes and items']);
                end
                L = L.data;
            elseif isnumeric(L)
                if ~(isvector(L) && length(L) == obj.numItems && min(L) >= 0 && ...
                     max(L) <= obj.numLabels && all(mod(L,1) == 0))
                    error(['Numeric label data must be a non-negative integer vector '
                           'and must match with respect to number of classes and items']);
                end
                L = L(:);
            end

            % Iterate over number of predictions
            Y = zeros(1,obj.numPredictions);
            if nargout >= 2
                A = zeros(obj.numLabels, obj.numLabels, obj.numPredictions);
            end
            for k = 1:obj.numPredictions
                I = obj.vectorFromIndexes(k);
                % Use only items with reference label and prediction data
                M = any(obj.data(:,I) > 0, 2) & (L > 0);
                % Multi-class AUC is computed as average of pair-wise AUC's
                Ak = MSAUCMatrix(obj.data(M,I), L(M));
                Y(k) = sum(Ak(:))/(obj.numLabels*(obj.numLabels-1));
                if nargout >= 2
                    A(:,:,k) = Ak;
                end
            end
        end
        
        function assert (obj)
            % validates score data (for each prediction the sub-row must be
            % either the zero vector or a vector of positive entries which
            % sum up to one.
            
            % Call superclass assert
            assert@MSPredictionData(obj);
            
            %verify that #columns is multiple of #classes
            if mod(obj.dataLength, obj.numClasses) ~= 0
                error(['The number of columns of the score data must be'...
                    ' a multiple of the number of classes'])
            end
            %verify that entries are in [0..1]
            if any(obj.data(:) < 0) || any(obj.data(:) > 1)
                error('The data entries must be in 0..1')
            end
%             %verify that for every score that the sum of rows is either 0
%             %or 1
%             for i=1:obj.numPredictions
%                 scores=obj.getPredictionData(i);
%                 rowSum=sum(scores,2);
%                 if ~all(abs(rowSum-0)<10^(-6) | abs(rowSum-1)<10^(-6))
%                     error(['The data does not represent a discrete pro'...
%                         'bability distribution'])
%                 end
%             end
        end   
    end
    
    methods(Access=protected)        
        function outData=aggregateMean(obj, scoreData, S, mask, numLabels)
            % aggregate scores according to mean criterium
            outData=zeros(numLabels, obj.numClasses);
            for i=1:obj.numClasses
                % Compute data per column
                maxLabelUsed=max(S(mask));
                outData(1:maxLabelUsed,i) = accumarray(S(mask), scoreData(mask,i),[],@mean);
            end
        end
        
        function outData=aggregateMedian(obj, scoreData, S, mask, numLabels)
            % aggregate scores according to median criterium
            outData=zeros(numLabels, obj.numClasses);
            for i=1:obj.numClasses
                % Compute data per column
                maxLabelUsed=max(S(mask));
                outData(1:maxLabelUsed,i) = accumarray(S(mask), scoreData(mask,i),[],@median);
            end
            % normalize discrete probability (as it is not invariant to median aggregation)
            obj.normalizeData(outData);
            
        end

        function indexVector=vectorFromIndexes(obj,index)      
            indexVector=zeros(1,length(index)*obj.numClasses);
            for i=1:length(index)
                indexVector((i-1)*obj.numClasses+1:...
                            (i-1)*obj.numClasses+obj.numClasses)=...
                            obj.numClasses*(index(i)-1)+1:...
                            obj.numClasses*index(i);
            end
        end
        
        function scoresData=select_impl(obj,index)
            % gets the score data for specified index 
            
            % gets the data (index validity is already checked)
            data=obj.getPredictionData(index, false);
            scoresData=MSScoreData(obj.labels,data,obj);
        end
        
        function value=getNumPredictions(obj)
            % gets number of predictions
            value=obj.dataLength/obj.numClasses;
        end
        
        function D = newFromTemplate (obj, data)
            % Create a new MSLabelData object with specified data
            D = MSScoreData(obj.labels, data);
        end
        
        function updateLabels_impl (obj, fwdLabelMap, bwdLabelMap)
            % Apply mapping to data matrix, unless mapping vectors equal (1:n)
            if any(fwdLabelMap ~= (1:length(fwdLabelMap))') || ...
               any(bwdLabelMap ~= (1:length(bwdLabelMap))')
                nc1 = length(fwdLabelMap);
                nc2 = length(bwdLabelMap);
                np = obj.numPredictions;
                newData = zeros(obj.numItems, nc2*np, 'like', obj.data);
                for k = 1:nc2
                    if bwdLabelMap(k)
                        i1 = (0:np-1)*nc1+bwdLabelMap(k);
                        i2 = (0:np-1)*nc2+k;
                        newData(:,i2) = obj.data(:,i1);
                    end
                end
                obj.data = newData;
            end
        end
    end
    
    methods(Static,Access=private)    
        function normalizeData(data)
            nrows=size(data,1);
            for j=1:nrows                
                try
                rowSum=sum(data(j,:));
                catch e
                    e.message
                end
                if rowSum~=0 && rowSum~=1
                    data(j,:)=data(j,:)/rowSum;
                end            
            end
        end
        function labelVector=score2LabelMax(data)
            % computes maximum per column and assigned the index
            [value,index]=max(data,[],2);
            % if maximum is 0 then leave spot without class asignment
            labelVector=index.*sign(value);
        end
        function labelVector=score2LabelProb(data)
            % computes an index according to probability distribution
            rows=size(data,1);
            labelVector=zeros(rows,1);
            for i=1:rows
                if sum(data(i,:))~=0
                    labelVector(i)=find(histc(rand,cumsum([0 data(i,:)])));
                end
            end
        end
    end
end