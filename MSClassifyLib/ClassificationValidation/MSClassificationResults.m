classdef MSClassificationResults < matlab.mixin.Copyable
 % This class stores the results for a classification task. A
 % classification task consists on performing the classification pipeline
 % with different preprocessing sequences (normalization, mean filtering,
 % spacial denoising, ...), feature extraction methods (Roc, PCA, NMF, ...),
 % map modifications (basis sorting, basis baseline removal, ...), and
 % classification models (LDA, SVM, ...)
 %
 % Properties:
 %  Properties consistent in cell array of structs containg data generated 
 %  using the whole dataset (field <full>) or specific data for cv-segment 
 %  (field <cv>):
 %  * predictionResults: predictions (based either in label or score data)
 %  * featureData: feature data
 %  * maps: feature maps
 %  * mapsModified: modifications of feature maps
 %  * classifModels: classification models
 %
 %  cvSegments:  cell array of MSLabelData objects containing the
 %      partition of data in CV segments 
 %  dataPartition: MSDataPartition used for the evaluation of the
 %      classifcation performance. Contains the ground truth in the property
 %      <classes>. Properties 'individuals', 'samples' and 'sampleSets'
 %      must contain the corresponding labels if methods 'aggregate' and
 %      'ROC4Aggregation' are called with a string input as aggregation
 %      reference.
 %  preprocessings: cell array with list of MSPreprocessingTrigger objects
 %      used to replicate the preprocessing pipeline using new data.
 %  numFeatures: cell array with list of number of features used during
 %      classification.
 %  errorReport: string containing errors during the classification study
 %  nErrors: count errors occurring during the classification study
 %  predictionType: string identifying the type of prediction. It can be
 %      one of the three values: 'label', 'scores', 'none'
 %  resultsTable: table object where the first columns are selectors,
 %      identifiers for the type of pipeline step and the rest of the columns
 %      contain the indexes where the objects corresponding to each column 
 %      are stored in the corresponding array property.
 %  nSelectors: Number of row selectors
 %  mainProcessSelectors: name of process selectors used for ordering the
 %      table columns after an append operation (just for reading convenience)
 %  mainStorageIndexes: name of storage columns
 %  numAppendedNewRegions: In case the method .appendNewRegions has been
 %      invoked this is not empty and shows how many regions have been
 %      appended to the original one for the given selector combination
 %  appendedDPartitionLookupTable: In case the method -appendNewRegions has
 %      been invoked this is not empty and shows which DPartition objects
 %      are actually appended versions of existing DPartition objects, e.g.
 %      if DPartition object 4 is an appended version of DPartition objects
 %      [1 2 3] in that order, this very property would save this
 %      information, i.e. in the 4th cell the vector [1 2 3] would be
 %      accessible
 %      --Remark: This implementation is not memory efficient. One could
 %      specify a get.Method for DPartition objects and execute the
 %      appending process on the fly with the information available in this
 %      very property. However, usual user cases would result in thousands
 %      of calls when calling .performance / .performanceRanking methods
 %      creating a lot overhead. On the other hand it is expected that the
 %      number of appended new regions, or their different combinations, is
 %      fairly small, such that the amount of redudantly saved combined
 %      DPartition objects is not a crucial factor.
 %      For now this property serves as user information and avoids
 %      multiple occurences or multiple concatenations of appended
 %      DPartition objects within single calls of .appendNewRegions
 %      by checking if a certain concatenation already
 %      exists before appending.
 %  
 % Properties (Dependent)
 %  nCombinations: number of different pipeline steps combinations
 %  selectorNames: name of selector columns
 %  storageNames: name of storage columns
 %
 % Methods
 %  apply: Classifies new data for the classification pipeline combinations
 %         given as input
 %  append: Appends one MSClassificationResults object to the other
 %  addSelector: Adds a new selector column (Useful to differentiate two
 %      previously appended objects
 %  select: Returns the indexes of given selector patterns 
 %  score2Label: Converts a score-based result into a label-based result
 %  prediction: Returns the prediction labels corresponding to a 
 %      given combination of processes
 %  findObject: Returns the object corresponding to a given combination of
 %      processes
 %  performanceRanking: Obtains the list of process combinations which rank
 %      highest/lowest according to a given performance criterium
 %  performance: Returns a struct with accuracy, kappa, sensitivity and 
 %      specificity per class corresponding to a given combination of
 %      processes
 %  plotPerformance: Graphics the accuracy, kappa, sensitivity and 
 %      specificity per class corresponding to a given performance struct
 %  aggregate: Returns a new MSClassificationResults corresponding to a
 %      certain aggregation level (e.g., patient level, core level, etc.)
 %  ROC4Aggregation: Obtains a list of MSROCCurve objects corresponding to
 %      a given aggregation level (e.g., patient level, core level, etc.) 
 %      and a specific combination of processes.
 %  predictionCombinations: Returns a list of valid combinations of processes
 %      according to the input selector pattern 
 %  printMeasureList: Print to the standard output a given list of
 %      combinations and associated performance measures
 %
 % MSClassificationResults uses the handle semantic, i.e. when 
 % assigning an object of this class to a variable, only a reference to
 % the original object is copied. Use the copy method to create a deep copy.
    properties (SetAccess = {?MSClassificationTask})        
        % cell arrays and arrays:
        predictionResults;
        preprocessings
        maps;        
        mapsModified; 
        featureData;        
        classifModels;        
        dataPartition;
        numFeatures;
        cvSegments;  
        errorReport = '';
        nErrors = 0;
        % information about the object structure
        predictionType;
        resultsTable;
        nSelectors;  
        mainProcessSelectors;
        mainStorageIndexes;
        numAdditionallyAppendedNewRegions = [];
        appendedDPartitionLookupTable = [];
    end  
    
    properties (Dependent)
        nCombinations;
        selectorNames
        storageNames
    end
   
    methods   % <dependent properties>
        function value = get.nCombinations(obj)
            value = size(obj.resultsTable,1);
        end        
        function value = get.selectorNames(obj)
            value=obj.resultsTable.Properties.VariableNames(1:obj.nSelectors);
        end
        function value = get.storageNames(obj)
            value=obj.resultsTable.Properties.VariableNames(obj.nSelectors+1:end);
        end
    end
    
    methods (Access = {?MSClassificationTask})  
        
        function obj = MSClassificationResults(string)
            % Constructor. It can be only invoked inside class MSClassificationTask
            % obj = MSClassificationResults
            % INPUT
            %   string: (optional) The string is displayed
            narginchk(0,1)
            if nargin > 0
                disp(string)
            end
        end
    end
    
    methods
        
        function results = apply(obj, newData, indexes, varargin)
            % Classifies new data according to the classification pipeline
            % indicated by the input process selector 'indexes'.
            % results = obj.apply(newData, indexes)
            % INPUT
            %   newData: MSMaldiData object containing the new data to be
            %            classified.
            %   indexes: array of indexes of combinations in resultsTable
            %            or a cell array of strings containing selector 
            %            patterns. If empty or not specified the
            %            classification is performed for each possible
            %            combination.
            %   varargin: name-value pair parameters:
            %       * predictionType: string indicating which kind of
            %           prediction should be computed (either 'scores' or
            %           'label').
            %       * saveFeatureData: logical value indicating whether the
            %           generated feature data during the classification
            %           process should be stored.
            %       * dPartition: MSDataPartitionObject containing among
            %           other labels the groundtruth used for evaluating
            %           performances.
            %       * resample: logical value indicating whether the data
            %           should be resampled in case of conflict between the 
            %           data mz resolution and the maps mz resolution. If
            %           resampling is set to false then in case of conflict
            %           an error will be thrown
            %       * cvIndex: Select CV model to apply
            %           cvIndex<0: related to the whole data 
            %           cvIndex==0: related to all cv-segments
            %           cvIndex>0: related to a specific cv-segment
            %           cvIndex=='cv': equivalent to cvIndex==0 (default)
            %           cvIndex=='full': equivalent to cvIndex==-1
            % OUTPUT
            %   results: MSClassificationResults object containing the
            %           results of the classification.           
            
            %check input
            narginchk(2,inf)
            if nargin < 3 || isempty(indexes)
                indexes = 1:obj.nCombinations;
            end  
          
            params=inputParser;
            params.addParameter('predictionType', 'label', @ischar);
            params.addParameter('saveFeatureData', false, @islogical);
            params.addParameter('dPartition',[],@obj.validateDataPartition);
            params.addParameter('resample',true, @islogical);
            params.addParameter('cvIndex','full');
            params.parse(varargin{:});
            predType = params.Results.predictionType;
            saveFD = params.Results.saveFeatureData;   
            dPartition = params.Results.dPartition;
            resample = params.Results.resample;
            cv = obj.validateCVIndex(params.Results.cvIndex);
            if cv==0
                error('You must specify a cv partition or the full data (-1 or full)')
            end
            % init output
            results = MSClassificationResults;  
            % get combinations
            indexes = obj.validateIndexes(indexes);
            % sort combinations according to order in resultsTable;
            indexes = sort(indexes);
            prep = [];
            mps = [];
            models = [];
            nFeatures = [];
            % get maps
            if ~isempty(obj.mapsModified)
                [mps, mapsIndex] = obj.findObject(indexes, 'map_modified',cv);
                if isempty(mps)
                    [mps, mapsIndex] = obj.findObject(indexes, 'map',cv);
                else
                    mask = mapsIndex > 0;
                    if any(~mask)
                        maxIndex = max(mapsIndex(mask));
                        [mps1, maps1Index] = obj.findObject(indexes(~mask), 'map', cv);
                        if ~isempty(mps1)
                            mapsIndex(~mask) = maps1Index + maxIndex*sign(maps1Index);
                            mps = horzcat(mps, mps1);
                        end
                    end 
                end
            elseif ~isempty(obj.maps)
                [mps, mapsIndex] = obj.findObject(indexes, 'map',cv);
            end
            if isempty(mps)
                error('No maps were stored for the input combinations')
            end
            % get models
            if ~isempty(obj.classifModels)
                [models, modelsIndex] = obj.findObject(indexes, 'model', cv);
            end
            if isempty(models)
                error('No models were stored for the input combinations')
            end
            mask = mapsIndex > 0 & modelsIndex > 0;
            if any(~mask)
                warning(['It is not possible to obtain a classification for '...
                        'some of the input combinations, as either maps or '...
                        'models were not stored during the training phase'])
            end
            % iterate through combinations....
            
            % masked indexes for which classification will be conducted
            indexes = indexes(mask);
            mapsIndex = mapsIndex(mask);
            modelsIndex = modelsIndex (mask);
            nComb = length(indexes);
            % init data needed to store generated objects
            lastMap = -1;
            predictions = cell(1,nComb);
            predictionIndex = zeros(nComb, 1);
            predictionLastIndex = 0;
            mData = newData;
            
            % get preprocessings
            if ~isempty(obj.preprocessings)
                [prep, prepIndex] = obj.findObject(indexes,'preprocessing');
                lastPrep = -1;                 
            end
            % get number of features arrays
            if ~isempty(obj.numFeatures)
                [nFeatures, nFeatIndex] = obj.findObject(indexes, 'number_features');
                results.numFeatures = nFeatures;
            end
            % initialize feature data array and indexes
            if saveFD
                FD = cell(1,nComb);
                FDindex = zeros(nComb, 1);
                FDLastindex = 0;
            end
            % start iteration....
            for i=1:nComb
                % preprocessing of the data
               if ~isempty(prep) && prepIndex(i)~=lastPrep
                   mData = newData.copy;
                   listPrepCurrent=prep{prepIndex(i)};
                   for j=1:length(listPrepCurrent)
                     listPrepCurrent(j).apply(mData);
                   end
                   lastMap = -1;
                   lastPrep = prepIndex(i);
               end
               % generate feature data
               currentMap = mps{mapsIndex(i)};
               % Can be the case when the case when cv > 0
               if iscell(currentMap) && length(currentMap)==1
                   currentMap=currentMap{1};
               end
               if mapsIndex(i)~=lastMap
                   lastMap = mapsIndex(i);
                   if  length(currentMap.mzVector) ~= mData.dataLength...
                      || norm(mData.mzVector - currentMap.mzVector) > 10^(-6)
                      if resample
                          mData.resample(currentMap.mzVector)
                      else
                          error(['The resolution of the input maldi data is '...
                              'different from the resolution of the stored map'...
                              ' and the results for this combinations will not'...
                              ' be included'])
                      end
                   end
                   featData = currentMap.map(mData);
                   
                   if saveFD
                       FDLastindex = FDLastindex + 1;
                       FD{FDLastindex} = featData;
                   end
               end
               if saveFD
                   FDindex(i) = FDLastindex;
               end
               % classify...
               
               % get current model(s)
               model = models{modelsIndex(i)};
               nNumFeatures = numel(model);
               classes = model{1}.classes;
               numClasses = numel(classes);
               
               if strcmp(predType,'scores')                             
                   prediction = ...
                         MSScoreData(classes,...
                                zeros(mData.numItems,...
                         nNumFeatures*numClasses), mData);
               else
                   prediction = ...
                        MSLabelData(classes,...
                                 zeros(mData.numItems,...
                                nNumFeatures), mData);
               end
               for k=1:nNumFeatures
                   if strcmp(predType,'scores')
                       [~, scores] =  model{k}.classify(featData); 
                       if ~isempty(scores)
                           prediction.data(:,numClasses*(k-1)...
                                            +1:numClasses*k) = scores.data;
                       else
                           prediction = [];
                           continue
                       end
                   else                      
                       kPred = model{k}.classify(featData);
                       prediction.data(:,k)=kPred.data;
                   end
               end
               if ~isempty(prediction)
                   predictionLastIndex = predictionLastIndex +1;
                   predictions{predictionLastIndex}.full = prediction;
                   predictionIndex(i) = predictionLastIndex;
               end
            end
            % Build results
            results.predictionResults = predictions(1:predictionLastIndex);
            % Build table....
            
            % include predictions
            t = array2table(predictionIndex);
            t.Properties.VariableNames = {'prediction'};
            % include number of features
            if ~isempty(nFeatures)
                t1 = array2table(nFeatIndex');
                t1.Properties.VariableNames = {'number_features'};
                t = [t t1];
            end
            % include feature data
            if saveFD
                FD = FD(1:FDLastindex);
                results.featureData = FD;
                t1 = array2table(FDindex);
                t1.Properties.VariableNames = {'feature_data'};
                t = [t t1];
            end
            % ensemble selectors and storage
            t = [obj.resultsTable(indexes,1:obj.nSelectors) t];
            results.resultsTable = t; 
            results.nSelectors = obj.nSelectors;
            results.predictionType = predType;
            if ~isempty(dPartition)
                if dPartition.numItems ~= newData.numItems
                    warning(['The input data partition cannot be stored as '...
                        'it does not contain the same number of items as '...
                        'the input maldi data'])
                else
                    results.setDataPartition(dPartition);
                end
            end
        end

        function append(obj,obj1) 
            %Appends one MSClassificationResults object to the other
            %obj.append(obj1)
            %   Input
            %   obj1: MSClassificationResults object to be appended to obj
            if ~isa(obj1,'MSClassificationResults')
                error('Input argument must be an MSClassificationResults object')
            end
            if (strcmp(obj.predictionType, 'scores') && ...
                 strcmp(obj1.predictionType,'label')) ||...
                 (strcmp(obj1.predictionType, 'scores') &&...
                 strcmp(obj.predictionType,'label'))
                disp(['The results must be based on the same kind of pre'...
                    'diction to be appended. Use the score2Label method'...
                    ' first in order two append label-based prediction results'])
                return             
            end
            % store two variables that will be modified in obj1 to keep
            % it unaltered after the process
            obj1Table = obj1.resultsTable;
            obj1nSelectors = obj1.nSelectors;
            % add missing selectors in both tables
            obj.addSelectorsFrom(obj1);
            obj1.addSelectorsFrom(obj);
            
            % get tables
            t=obj.resultsTable;
            t1=obj1.resultsTable;            
            % get number of selectors
            nS1=obj1.nSelectors;
            nS=obj.nSelectors;            
            % get selector column names
            varNames=t(:,nS+1:end).Properties.VariableNames;
            varNames1=obj1.resultsTable(:,nS1+1:end).Properties.VariableNames;
            
            % add columns in t1 from t. If columns exist in both tables
            % update the indexes and the corresponding property of obj
            for i=1:width(t)-nS
                index=find(strcmp(varNames{i},varNames1),1);
                if isempty(index) % the column does not exist in t1
                    obj1.addEmptyStorage(varNames{i});
                else % the none-zero entries of t1 must be incremented
                    mask=t1{:,index+nS1}>0;
                    maxIndex=max(t{:,i+nS});
                    obj1.resultsTable{mask,index+nS1}=obj1.resultsTable{mask,index+nS1}+maxIndex;
                    [~,propertyName] = obj.listFromColumn(varNames{i});
                    list = obj1.listFromColumn(varNames{i});
                    obj.(propertyName) = horzcat(obj.(propertyName),list);
                end
            end            
            % Finally add columns in t from t1
            for i=1:width(t1)-nS1
                index=find(strcmp(varNames1{i},varNames),1);
                if isempty(index) %only necessary to consider this case
                    obj.addEmptyStorage(varNames1{i});
                    [list,propertyName] = obj1.listFromColumn(varNames1{i});
                    obj.(propertyName) = list;
                end
            end
            
            % append table
            obj.resultsTable = vertcat(obj.resultsTable, obj1.resultsTable);
            % Leave obj1 unaltered
            obj1.resultsTable = obj1Table;
            obj1.nSelectors = obj1nSelectors;
        end
        
        function appendNewRegions(obj,obj1) 
            % Appends one MSClassificationResults object to the other
            % obj.appendNewRegions(obj1)
            %
            % The difference to .append is that no new rows or prediction
            % results are being appended as a whole, 
            % instead prediction results are being merged 
            % and possibly existing position informations discarded. 
            % Therefore existing selector names of both objects
            % are required to match, i.e. each row in obj1 has to exist in
            % obj exactly one time.
            %
            % This functionality is especially useful for sequentially 
            % appending validation sets with trivial groundtruth 
            % (e.g. only one label occurs in the data) in order to get an
            % AUROC score performance across all appended sets for specific
            % selector combinations.
            %
            % Remarks: it is not checked if numFeatures match, it might be
            %   user intended to append even without such a match,
            %   dataPartition objects are concatenated accordingly, and
            %   positions deleted there as well
            %
            %   Input
            %   obj1: MSClassificationResults object to be appended to obj
            if obj==obj1
              error(['It is prohibited to concatenate the object to itself '...
                'since handle object behavior would result in multiple '...
                'internal overwritings of certain properties']);
              % If the user really wants this strange behavior he can
              % create a hard copy of the results object with .copy
            end
            if ~isa(obj1,'MSClassificationResults')
                error('Input argument must be an MSClassificationResults object')
            end
            if (strcmp(obj.predictionType, 'scores') && ...
                 strcmp(obj1.predictionType,'label')) ||...
                 (strcmp(obj1.predictionType, 'scores') &&...
                 strcmp(obj.predictionType,'label'))
                disp(['The results must be based on the same kind of pre'...
                    'diction to be appended. Use the score2Label method'...
                    ' first in order two append label-based prediction results'])
                return             
            end
            if ~(obj.nSelectors == obj1.nSelectors && ...
                all(strcmp(obj.selectorNames,obj1.selectorNames)))
              error('Result object selector names do not match');
            end
            % get number of selectors
            nS=obj.nSelectors;                          
            
            % get tables
            t1=obj1.resultsTable(:,1:nS);
            t1c = table2cell( t1 );
            
%            t= obj.resultsTable(:,1:nS);
%            tc = table2cell( t );
%             rowCompareFun = @(x,x1)(all(strcmp(x,x1)));
%             wrappedCellArrayT = cell( size(t,1), 1 );
%             for i = 1:size(t,1)
%               wrappedCellArrayT{i} = tc(i,:);
%             end
            
            if isempty(obj.numAdditionallyAppendedNewRegions)
              obj.numAdditionallyAppendedNewRegions = zeros( obj.nCombinations, 1 );
            end

            for j = 1:size(t1,1)
              x1 = t1c(j,:);
              %hits = cellfun( @(x)rowCompareFun(x,x1), wrappedCellArrayT );
              hits = obj.select( x1{:} );
              %nHits = sum(hits);
              nHits = numel( hits );
              if nHits > 1
                error('Some selector name combinations in the object calling this method are not unique');
                % If this is a desired use case the behavior can be changed
              elseif nHits < 1
                error('Some selector name combinations of the object to be appended could not be found in the object calling this method');
              end
              %
              obj.numAdditionallyAppendedNewRegions(hits) = obj.numAdditionallyAppendedNewRegions(hits)+1;
              oldPrediction = obj.predictionResults{hits};
              newPrediction = obj1.predictionResults{j};
              obj.predictionResults{hits} = mergePredictions(oldPrediction,newPrediction);
              %
              
              % Concatenate data partition objects ( groundtruth ) if
              % necessary
              if j == 1
                numOldDPartitions = numel(obj.dataPartition);
                numNewDPartitions = numel(obj1.dataPartition);
                obj.dataPartition = [ obj.dataPartition obj1.dataPartition ];
                if isempty(obj.appendedDPartitionLookupTable)
                  obj.appendedDPartitionLookupTable = cell(1,numOldDPartitions);
                  for k = 1:numOldDPartitions
                    obj.appendedDPartitionLookupTable{ k } = k;
                  end
                end
                appendixDP = cell(1,numNewDPartitions);
                for k = 1:numNewDPartitions
                  appendixDP{k}= numOldDPartitions + k;
                end
                obj.appendedDPartitionLookupTable = [...
                  obj.appendedDPartitionLookupTable appendixDP];
              end
              oldDPartitionIdx = obj.resultsTable(hits,:).data_partition;
              newDPartitionIdx = obj1.resultsTable(j,:).data_partition;
              newDPartitionIdx = newDPartitionIdx + numOldDPartitions;
              currComb = [ oldDPartitionIdx newDPartitionIdx ];
              checkCurrComb = cellfun( @(x)(all(x==currComb)), ...
                obj.appendedDPartitionLookupTable, 'ErrorHandler',@(x,y)(false) );
              switch sum(checkCurrComb)
                case 0
                  % New entry required - concatenate DPartition objects
                  oldDPartition = obj.dataPartition{ oldDPartitionIdx };
                  newDPartition = obj.dataPartition{ newDPartitionIdx };
                  combinedDPartition = concatDPartitions( oldDPartition, newDPartition );
                  currDPartitionIdx = numel( obj.dataPartition ) + 1;
                  obj.dataPartition{currDPartitionIdx} = combinedDPartition;
                  obj.appendedDPartitionLookupTable{ currDPartitionIdx } = currComb;
                case 1
                  % Entry existing - no new obj.dataPartition entry
                  % necessary - should be the usual case!!
                  currDPartitionIdx = find( checkCurrComb );
                otherwise
                  error('Multiple hits for particular dataPartition combination according to lookup table - this should not happen');
              end
              %
              obj.resultsTable.data_partition(hits) = currDPartitionIdx;
              % Remark: a similar concatenation / update behavior could be
              % defined for the nFeatures storage - however, it is expected
              % that the user combines identical numbers of features or
              % knows what he / she is doing if the nFeatures of the
              % appended results do not match.
              % The prediction storage idx should not change, however,
              % technically the prediction it points to is changed.
            end
            obj.cleanUpUnusedDPartitionObjects;
            
        end
        
        function cleanUpUnusedDPartitionObjects(obj)
            totalDPObj = 1:numel(obj.dataPartition);
            usedDPObj = unique( obj.resultsTable.data_partition );
            unusedDPObj = setdiff( totalDPObj, usedDPObj );
            obj.dataPartition( unusedDPObj ) = cell(size(unusedDPObj));
        end
        
        function addSelector(obj, columnName, columnValue)
            % Adds a new selector in order to retrieve rows from the
            % combination table resultsTable (which is modified with an 
            % extra column)
            % obj.addSelector(columnName, columnValue)
            % INPUT
            %   columnName: Name of the selector (column in resultsTable)
            %   columnValue: cellarray of strings identifying the selector
            %       value for the different process combinations. It can 
            %       also be a string, in which case all the combinations
            %       have the same selector value (the input string).
            
            %input validation
            narginchk(2,3)
            if ~ischar(columnName)
                error('Input <columnName> must be a string')
            end
            if any(strcmp(columnName, obj.resultsTable.Properties.VariableNames))
                error(['There exists already a selector with name <'...
                     columnName '>. Try with another name'])
            end
            if nargin < 3
                columnValue='none';
            end
            % if columnValue is a string a column is added with the value
            % repeated for every combination
            if ischar(columnValue)
                c = cell(obj.nCombinations,1);
                c(:) = {columnValue};
                columnValue = c;
            else % in case of a cellarray check size of input
                if ~(iscellstr(columnValue) && norm(size(columnValue)-...
                                                [obj.nCombinations, 1])==0)
                    error(['Input <columnValue> must be a single-column '...
                        'cell array of %d elements'], obj.nCombinations)                        
                end
            end
            % convert column to table
            t = cell2table(columnValue);
            % add column name
            t.Properties.VariableNames = {columnName};
            % insert column as last selector
            obj.resultsTable=[obj.resultsTable(:,1:obj.nSelectors) t...
                                  obj.resultsTable(:,obj.nSelectors+1:end)];
            %increase nSelectors property
            obj.nSelectors = obj.nSelectors + 1;
        end                   

        function [object, matchingIndexes]=findObject(obj,indexes,...
                                                    objectType, cvIndex, once)
            % Returns an object corresponding to a given combination of
            % processes. The object can be any of those for which a storage
            % column in the resultsTable (not a selector column) can be
            % found.
            % [object, matchingIndexes] = obj.findObject(indexes, objectType, cvIndex, once)
            % INPUT
            %   indexes: array of indexes of combinations in resultsTable
            %    or a cell array of strings containing selector patterns
            %   objectType: string indicating which type of object should
            %    be returned, e.g. map, map_modified, model, etc. The string
            %    must be one of the resultsTable storage column names
            %   cvIndex: numerical scalar or string
            %       cvIndex<0: related to the whole data 
            %       cvIndex==0: related to all cv-segments
            %       cvIndex>0: related to a specific cv-segment
            %       cvIndex=='cv': equivalent to cvIndex==0 (default)
            %       cvIndex=='full': equivalent to cvIndex==-1
            %   once: string 'once' or logical value indicating whether
            %       only the first result from indexes will be taken into
            %       account. If <once> then the expected output will never 
            %       be a one-element cell array
            % OUTPUT
            %   object: cell array (or object in case 'once') of objects.
            %   matchingIndexes: index array with as many elements as
            %       requested combinations indicating in which position of
            %       the array object the related object is found
            
            
            %input validation
            narginchk(3,5)
            if nargin < 5
                once = false;
            elseif ischar(once) && strcmp(once,'once')
                once = true;
            else
                assertLogicalScalar(once, 'once');
            end
            if nargin<4||isempty(cvIndex)
                cvIndex=0;                
            else
                cvIndex=obj.validateCVIndex(cvIndex);
            end
            % compute table indexes
            indexes=obj.validateIndexes(indexes, once);
            if isempty(indexes)
                disp('No results for this combination')
                object=[]; matchingIndexes=[];
                return
            end
            % decide which list of objects is to be used
            list = obj.listFromColumn(objectType);
            
            if isempty(list)
                disp('The type of requested object was not stored for the input combinations')
                object=[]; matchingIndexes=[];
                return
            end
            % find indexes of objects in table
            objectIndexes=obj.resultsTable(indexes,:).(objectType);
            % mask of positive indexes
            mask = objectIndexes > 0;
            if ~any(mask)
                disp('The type of requested object was not stored for the input combinations')
                object=[]; matchingIndexes=[];
                return
            end
            matchingIndexes=zeros(1,length(objectIndexes));
            % compute unique coeficients
            [uniqueIndexes,~,matchingIndexes(mask)]=unique(objectIndexes(mask));
            % get list of objects
            object=list(uniqueIndexes);
            nObjects=numel(object);
            % filter object according to cvIndex, unless a preprocessing is requested 
            if ~(any(strcmp(objectType,{'preprocessing','data_partition','cv_segment','number_features'})))
                for i=1:nObjects
                    if cvIndex==-1
                        if isfield(object{i},'full')
                            object{i} = object{i}.full;
                            
                        else
                            object{i} = [];
                            matchingIndexes(matchingIndexes==i)=0;
                        end
                    else
                        if ~isfield(object{i},'cv')
                            object{i} = [];
                        else
                            %Includes the case where more than one
                            %cvSegment object is required
                            if cvIndex > 0 || length(cvIndex)>1 
                                try
                                    object{i} = object{i}.cv(:,cvIndex);
                                catch % the requested cvIndex exceeds the length of the cv array
                                    object{i} = [];
                                end
                            else
                                object{i} = object{i}.cv;
                            end
                        end
                    end
                end
                % Eliminate empty objects (caused by not having the required fields)
                nonEmptyObj=~cellfun(@isempty,object);
                if isempty(find(nonEmptyObj,1))
                    disp('The type of requested object was not stored for the input combinations')
                    object = []; matchingIndexes=[];
                    return
                else
                    object = object(nonEmptyObj);
                    matchingIndexes = obj.removeIndexes(matchingIndexes, nonEmptyObj);
                end
            end

            % if requested only one object extract it from the outer cell
            % array
            if once && iscell(object) && numel(object)==1
                object=object{1};                
            end
        end  
  
        function indexes = select(obj,varargin)
            % Return a list with the indexes of those combinations matching
            % the input selector's patterns
            % indexes = obj.select(patternSelector1, patternSelector2,...)
            % INPUT
            %   varargin: strings containing the selector patterns. The
            %   order of the patterns is according to the order of the
            %   selector columns in resultsTable
            % OUTPUT
            %   indexes: list with the indexes of the requested
            %     combinations in resultsTable
            
            % check input
            narginchk(1,obj.nSelectors+1);
            if nargin == 1 || (numel(varargin)==1 && isempty(varargin{1}))
                indexes = 1:obj.nCombinations;
                return
            elseif ~iscellstr(varargin)
                error('All input arguments must be strings')
            end
            % get column names
            colNames = obj.resultsTable.Properties.VariableNames;
            % matching per column and return intersection
            indexes = true(obj.nCombinations,1);
            for i=1:length(varargin)
                if ~isempty(varargin{i})
                    match = regexp(obj.resultsTable.(colNames{i}),varargin{i},'once');
                    indexes=indexes & ~cellfun('isempty',match);
                end
            end
            % obtain indexes from logical vector
            indexes=find(indexes);
        end                          
        
        function [predictions,indexes]=prediction(obj, indexes, cvIndex, once, joinCVSegments)
            % Returns the prediction labels corresponding to a specified
            % process combination.
            % pLabels=obj.predictionLabels('tic','nmf','none','SVM', true)
            % pLabels=obj.predictionLabels('tic',1,'none',2)
            % INPUT
            %   indexes: cell array with processes search patterns or
            %   numerical array with indexes of process in obj.resultsTable
            %   cvIndex: numerical scalar or string
            %       cvIndex<0: related to the whole data 
            %       cvIndex==0: related to all cv-segments
            %       cvIndex>0: related to a specific cv-segment
            %       cvIndex=='cv': equivalent to cvIndex==0 (default)
            %       cvIndex=='full': equivalent to cvIndex==-1
            %   once: string 'once' or logical value indicating whether
            %       only the first result from indexes will be taken into
            %       account. If <once> then the expected output will never 
            %       be a one-element cell array
            %   joinCVSegments: flag that decides if a joined output is
            %       generated (default) or if cv segment predictions will
            %       be returned separately (in a cell array)
            % OUTPUT
            %   predictions: prediction label
            %   indexes: optional output argument returning the indexes
            
            %input validation
            narginchk(2,5)
            if nargin < 5
                joinCVSegments = true;
            end
            if nargin < 4
                once = false;
            elseif ischar(once) && strcmp(once,'once')
                once = true;
            else
                assertLogicalScalar(once, 'once');
            end
            if nargin<3||isempty(cvIndex)
                cvIndex=0;                
            else
                cvIndex=obj.validateCVIndex(cvIndex);
            end
            % obtain (and validate) indexes in results table
            indexes=obj.validateIndexes(indexes, once);
            
            if isempty(indexes)
                disp('No combinations were generated for the specified indexes');
                predictions=[];
                return;
            end
            
            %get prediction indexes
            predIndexes = obj.resultsTable.prediction(indexes);
            %mask of valid indexes (>0)           
            indexMask=predIndexes>0;
            if ~any(indexMask)
                disp('No predictions were stored for this process combination');
                predictions=[];
                return;
            end
            % obtain cv_segment list
            if isscalar(cvIndex) && cvIndex > 0 || length(cvIndex) > 1
                [cvSegm, indexCVSegm] = obj.findObject(indexes, 'cv_segment');
                if isempty(cvSegm)||~all(indexCVSegm)
                    error(['A cv_segment is needed to compute the predictions but '...
                          'it was not stored for some of the combinations'])
                end
                indexCVSegm=indexCVSegm(indexMask);
            end
            % initialize predictions
            predictions=cell(length(indexes),1);
            % obtain struct with predictions
            predictionData=obj.predictionResults(predIndexes(indexMask));
            emptyObject=true;
            for i=1:sum(indexMask)            
                if all(cvIndex==-1) % return predictions for the whole dataset
                    if isfield(predictionData{i},'full')
                        predictionData{i}=predictionData{i}.full;
                        emptyObject=false;
                    end
                else
                    if ~isfield(predictionData{i},'cv')
                        continue
                    else
                        emptyObject=false;
                    end
                    if all(cvIndex==0) % returns a cell array with predictions
                                  % per cv-segment
                        predictionData{i}=predictionData{i}.cv;
                        
                    else % returns the prediction for the specified cross-validation segment   
                        % includes the case where more than one cv is
                        % required
                        predictionD = cell(1,length(cvIndex));
                        for j=1:length(cvIndex)
                            if ~iscell(predictionData{i}.cv)
                                mask=cvSegm{indexCVSegm(i)}.data~=cvIndex(j);
                                if isa(predictionData{i}.cv,'MSLabelData')
                                    predictionD{j}=MSLabelData(predictionData{i}.cv.labels...
                                        ,predictionData{i}.cv.data, predictionData{i}.cv);
                                    predictionData{i}.data(mask,:)=0;
                                else
                                    predictionD{j}=MSScoreData(predictionData{i}.cv.labels...
                                        ,predictionData{i}.cv.data, predictionData{i}.cv);
                                    predictionData{i}.data(mask,:)=0;
                                end
                            else
                                predictionD{j}=predictionData{i}.cv{cvIndex(j)};
                            end
                        end
                        
                        if length(cvIndex)==1
                            predictionData{i}=predictionD{1};
                        else
                            if joinCVSegments
                                %join all the predictions
                                predictionData{i} = MSPredictionData.join(predictionD);
                            else
                                predictionData{i} = predictionD;
                            end
                        end
                    end
                end
            end
            if emptyObject
                disp('No predictions were obtained for the specified combinations')
                predictions=[];
                return
            end
            predictions(indexMask)=predictionData;
            if numel(predictions)==1 && once
                predictions=predictions{1};
            end
        end
        
        function performanceResults = performance(obj, indexes, cvIndex, once, joinCVSegments)
            % Returns a structure with the performance of the classification
            % corresponding to a specified process combination.
            % performanceResults=obj.performance('tic','nmf','none','SVM', true, false)
            % INPUT
            %   indexes: cell array with processes search patterns or
            %   numerical array with indexes of process in obj.resultsTable
            %   cvIndex: numerical scalar or string
            %       cvIndex<0: classification using whole data for training 
            %                   testing
            %       cvIndex==0:combined performance for all cv-segments
            %       cvIndex>0: performance for the specified cvIndex
            %       cvIndex == 'all':
            %       cvIndex=='cv': equivalent to cvIndex==0 (default)
            %       cvIndex=='full': equivalent to cvIndex==-1
            %   once: string 'once' or logical value indicating whether
            %       only the first result from indexes will be taken into
            %       account. If <once> then the expected output will never 
            %       be a one-element cell array
            %   joinCVSegments: flag that decides if a joined output is
            %       generated (default) or if cv segment predictions will
            %       be returned separately (in a nested cell array)

            % OUTPUT
            %   performanceResults: cellarray of structure with the following fields:
            %      -accuracy: list of positive doubles
            %      -kappa: list of positive doubles
            %      -sensitivity: matrix of positive doubles
            %      -specificity: matrix of positive doubles
            %      -avgSensitivity: list of positive doubles
            %      -avgSpecificity: list of positive doubles

            %input validation
            narginchk(2,5)
            if nargin < 5
              joinCVSegments = true;
            end
            if nargin < 4
                once = false;
            elseif ischar(once) && strcmp(once,'once')
                once = true;
            else
                assertLogicalScalar(once, 'once');
            end
 
            if nargin < 3||isempty(cvIndex)
                cvIndex=0;                
            else
                cvIndex = obj.validateCVIndex(cvIndex);
            end
            % obtain predictions and suppress joining cv segments at this stage
            [prediction,indexes]=obj.prediction( indexes, cvIndex, false, false );
            % obtain ground truth
            [dPartition, indexPartition] = obj.findObject(indexes, 'data_partition');
            if isempty(prediction)
                %no results were generated for such combination
                disp('No predictions were generated for the given combination')
                performanceResults=[];
                return
            end
            % Obtain performance struct for each prediction
            nPred=numel(prediction);
            performanceResults=cell(nPred,1);
                         
            for k=1:nPred
              groundTruth = dPartition{ indexPartition(k) }.classes;
              if ~isempty( prediction{k} )
                if joinCVSegments
                  predictionKJoined = MSPredictionData.join( prediction{k} );
                  performanceResults{k} = obj.performanceFromLabel( predictionKJoined, groundTruth );
                else
                  nCV = numel( prediction{k} ); % Cell array with nCV cells
                  performanceResults{k} = cell(nCV,1);
                  if nCV == 1
                    performanceResults{k}{nCV} = obj.performanceFromLabel( prediction{k}, groundTruth );
                  else
                    for iCV = 1:nCV
                      performanceResults{k}{iCV} = obj.performanceFromLabel( prediction{k}{iCV}, groundTruth );
                    end
                  end
                end
              end
              % Handle cases cvSelf if relevant
              if isfield( obj.predictionResults{indexes(k)}, 'cvSelf' ) && ...
                    ~isempty( obj.predictionResults{indexes(k)}.cvSelf ) && ...
                    ( isa( obj.predictionResults{indexes(k)}.cvSelf, 'MSScoreData') || ...
                    (iscell( obj.predictionResults{indexes(k)}.cvSelf ) && ...
                    isa( obj.predictionResults{indexes(k)}.cvSelf{1}, 'MSScoreData')) )
                nCV = numel( obj.predictionResults{indexes(k)}.cvSelf ); % Cell array with nCV cells
                if nCV == 1
                  [labelPredictionkAfterThr, optimalThresholds] = ...
                    MSScore2LabelViaThreshold( prediction{k}, obj.predictionResults{indexes(k)}.cvSelf, groundTruth );
                  labelPredictionkAfterThrJoined = MSPredictionData.join( labelPredictionkAfterThr );
                  CVTrainedThrPerf = obj.performanceFromLabel( labelPredictionkAfterThrJoined, groundTruth );
                  performanceResults{k}.CVTrainedThrPerf = CVTrainedThrPerf;
                  performanceResults{k}.CVTrainedThrPerf.usedThresholds = optimalThresholds;
                else
                  for iCV = 1:nCV
                    [labelPredictionkAfterThr, optimalThresholds] = ...
                      MSScore2LabelViaThreshold( prediction{k}{iCV}, obj.predictionResults{indexes(k)}.cvSelf{iCV}, groundTruth );
                    CVTrainedThrPerf = obj.performanceFromLabel( labelPredictionkAfterThr{1}, groundTruth );
                    performanceResults{k}{iCV}.CVTrainedThrPerf = CVTrainedThrPerf;
                    performanceResults{k}{iCV}.CVTrainedThrPerf.usedThresholds = optimalThresholds;
                  end
                end
              end
            end
            
            if numel(performanceResults)==1 && once
                performanceResults=performanceResults{1};
            end
        end   
        
        function [ranking,  value]=performanceRanking(obj,funHandle,cvIndex,...
                                                indexes, numSelected)
        % Returns a list of combinations ordered according to a given
        % performance criterium.
        % obj.performanceRanking(@(x) mean(x.accuracy),-1, 10)
        % obj.performanceRanking(@(x) x.accuracy(end),false, -10,combinations)
        % INPUT
        %   funHandle: function handle which receives either:
        %      *a structure
        %       representing the performance of the classification and
        %       returns a numerical value used for ranking (higher=better).
        %       The performance structure has the following fields:
        %        -accuracy: list of positive doubles
        %        -kappa: list of positive doubles
        %        -sensitivity: matrix of positive doubles
        %        -specificity: matrix of positive doubles
        %        -avgSensitivity: list of positive doubles
        %        -avgSpecificity: list of positive doubles
        %        Function call:value=funHandle(performance_struct)
        %       *indexes indicating specific combination and an integer
        %        cvIndex to specify the type of predictions (full data,
        %        combination of cv-performances or individual 
        %        cv-performance are considered.
        %        Function call:value=funHandle(iPP,iFE,iMM,iC,allData)
        %   cvIndex: numerical scalar or string
        %       cvIndex<0: classification using whole data for training 
        %                   testing
        %       cvIndex==0:combined performance for all cv-segments
        %       cvIndex>0: performance for the specified cvIndex
        %       cvIndex=='cv': equivalent to cvIndex==0 (default)
        %       cvIndex=='full': equivalent to cvIndex==-1
        %   numSelected: (optional)integer indicating how many combinations
        %       should be included. If numSelected>0 the best combinations
        %       are shown, otherwise the worst combinations are shown. 
        %       (default=total number of valid combinations)
        %   combinations: (optional) a cell array where each element is a 
        %       cell array containing a specific combination. 
        %       (default=all valid combinations)
        % OUTPUT
        %   ranking: a cell array where each element is a cell array 
        %       containing a specific combination. The elements are ordered
        %       according to the specified performance criterium
        %   value:(optional) a list with the value for each combination
        %       used to do the ranking
        
            %input validation
            narginchk(2,5)
            if nargin < 4 ||isempty(indexes)
                % if combinations not specified, use all valid combinations
                indexes = obj.predictionCombinations;          
            else
                indexes = obj.predictionCombinations(indexes);
            end
            if isempty(indexes)
                ranking = [];
                value = [];
                return
            end
            nIndexes=length(indexes);
                       
            if nargin < 5||isempty(numSelected)
                % if numSelected not specified consider all combinations
                numSelected=nIndexes;
            elseif ~isscalar(numSelected)||~isnumeric(numSelected)||abs(numSelected)<1
                error('Input <numSelected> must be an integer')
            else
                numSelected=sign(numSelected)*min(floor(abs(numSelected)),nIndexes);
            end
            if nargin<3||isempty(cvIndex)
                cvIndex=0;                
            else
                cvIndex = obj.validateCVIndex(cvIndex);
            end         
            performanceResults=cell(nIndexes,1);
            try
                % obtain performances
                performanceStruct=obj.performance(indexes, cvIndex);
                if isempty(performanceStruct)
                    disp('No predictions were computed for the specified cvIndex')
                    ranking=[];value=[];
                    return
                else
                    for i=1:nIndexes
                        if ~isempty(performanceStruct{i})
                            performanceResults{i}=funHandle(performanceStruct{i});
                        end
                    end
                end
            catch
                try
                    for i=1:nIndexes
                        performanceResults{i}=funHandle(obj,indexes(i),cvIndex);
                    end
                catch ME
                    msg=ME.message;
                    error(['Something went wrong. Either the function '...
                       'handle <funHandle> or the list of combinations'...
                       ' <combinations> is wrongly defined. More infor'...
                       'mation about the error: ' msg])
                end
            end
            mask=~cellfun('isempty',performanceResults); %Use the mask
            numSelected=min(numSelected,sum(mask)); % necessary in case of empty results
            if nIndexes==1
                index=1;
            else
                %order of elements in ascending order
                try
                    [~,index]=sort(cell2mat(performanceResults),1);
                catch
                    error(['Different number of performance measures has '...
                        'been generated for different combinations. Try with'...
                        ' a different ranking function that generates the'...
                        ' same number of measures per combination'])
                end
                if numSelected<0
                    %return elements in ascending order
                    index=index(1:-numSelected);
                else
                    %return elements in descending order
                    index=index(end:-1:end-numSelected+1);
                end
            end
            maskedIndexes=indexes(mask);
            ranking=maskedIndexes(index);
            maskedResults=performanceResults(mask);
            value=maskedResults(index);
            %print results
            obj.printMeasureList(ranking, value);
        end          
                
        function combinations=predictionCombinations(obj, indexes)
            % Returns a cell array with indexes of all (or selected from 
            % specified indexes) combinations for which a
            % prediction was computed. 
            % combinations=obj.predictionCombinations(numeric,PP,FE,MM,C)
            % INPUT
            %   indexes: cell array with processes search patterns or
            %   numerical array with indexes of process in obj.resultsTable
            % OUTPUT
            %   combinations: array with indexes of the obj.resultsTable
            %       for which predictions are stored
            
            %input validation
            if ~strcmp(obj.storageNames,'prediction')
                disp('No predictions were stored for this task')
                combinations = [];
                return
            end
            narginchk(1,2)
            if nargin<2
                indexes = (1:obj.nCombinations);
            else
                indexes = obj.validateIndexes(indexes);
            end
            if isempty(indexes)
                disp('No predictions were stored for these combinations')
                combinations = [];
                return
            end
            combinations = indexes(obj.resultsTable.prediction(indexes)>0);          
        end            
        
        function rTable = printMeasureList(obj,indexes,results)
            % Prints to the standard output a list with the strings 
            % corresponding to a process combination and a result (measure)
            % associated to each combination
            % obj.printMeasureList(indexes,results)
            % INPUT
            %   indexes: cell array where each element is a cell array of
            %       integers or strings indicating a specific process
            %       combination
            %   results: cell array with length(results)=length(indexes) 
            %       where each element is a vector of measures associated
            %       with each combination
            if numel(obj.numFeatures) == 1 % numFeature-Vector unique
              rTable=[obj.resultsTable(indexes,1:obj.nSelectors)...
                      array2table(cell2mat(results))];
            else
              rTable=[obj.resultsTable(indexes,1:obj.nSelectors)...
                      array2table(cell2mat(results))];
            end
            disp(rTable);
        end                        
        
        function outResults=aggregate(obj, S, aggregation)
          % Create new MSClassificationResults object
          %   representing label values aggregated over all items with the same
          %   label value in segment label object S (MSLabelData).
          %
          %   By default, aggregation is done by selecting the positive label
          %   value that occurs most often among all values with the same S
          %   label value ('max' aggregation). In case more than one label
          %   value occurs the maximum number of times, the lowest value is
          %   selected.
          %
          % outResults = obj.aggregate(S, aggregation): Use specified
          %   aggregation method to aggregate labels over segments. aggregation
          %   may be one of the following:
          % INPUT
          %
          %  
          % ---------------------------------------------------------------
          %  aggregation:
          %     A (1,2)-cell with each entry being one of the following
          %
          %   If Labels:
          %     'max' - Select positive label value occuring most often within
          %             each segment. If not unique, select the lowest such
          %             label. Items with label value == 0 are ignored.
          %     'unique' - Select 0 if all label values within segment equal 0,
          %                select NaN if >= 2 different values > 0 occur within segment
          %                  (this raises a warning),
          %                otherwise select the unique positive value occuring in segment.
          %  If Scores:
          %     'median' or 'mean', compare MSScoreData
          % 
          %  For both PredictionDataTypes:
          %     fn - A function handle fn that is called for each segment with
          %          two arguments. The first argument is a vector containing
          %          the number of occurences of all possible positive label 
          %          values, the second argument is the total number of items
          %          in the segment (including items with label value == 0).
          %          fn must return the selected label value as a scalar.
          %
          %  The first cell entry indicates the intended aggregation for
          %  the prediction data, the second cell entry indicates the
          %  intended aggregation for the groundtruth data - usually the
          %  latter is MSLabelData even when the first is MSScoreData and
          %  usually the default, i.e. aggregation{2} == [], selecting
          %  'max' is reasonable. To further enhance compatibility with old
          %  code, if the input for aggregation is (only one) function
          %  handle or of type char, it will be transformed to a (1,2)-
          %  cell with the 2nd entry being empty.
          %
          % ---------------------------------------------------------------
          %
          %   S: MSLabelData object indicating aggregation level
          %  OUTPUT
          %     outResults: MSClassificationResults object representing
          %         results for aggregation
            
            %input validation
            narginchk(1,3)
            if nargin < 2 ||  isempty(S)
              S = []; % uses trivialAggregation defined in MSPredictionData
            end
            if nargin<3
                aggregation = cell(1,2);
            elseif numel( aggregation ) == 1 || ischar( aggregation )
              aggregationC = cell(1,2);
              aggregationC{1} = aggregation;
              aggregation = aggregationC;
            end
            
            outResults = obj.copy;
            if ~isempty(obj.predictionResults)                
                indexes = obj.predictionCombinations;
                nComb = length(indexes);  
                if ischar(S)
                    [dPartition, indexPartition] = obj.findObject(indexes, 'data_partition');
                else
                    aggrLabel = S;
                end
                % mask of predictions which are not possible to be computed
                % because the aggregation label is not specified for the
                % given combinations
                predictionMask = false(1, nComb);
                pColumn = obj.resultsTable.prediction;
                for k=1:nComb                                       
                    % if trained with all data aggregate the
                    % prediction labels corresponding to that
                    % training
                    pResults=obj.predictionResults{pColumn(indexes(k))};                    
                    if ischar(S)
                        if ~indexPartition(k)
                            predictionMask(k) = true;
                            continue
                        end
                        aggrLabel=dPartition{indexPartition(k)};
                        
                        if ~any(strcmp(S,fieldnames(aggrLabel)))...
                                ||isempty(aggrLabel.(S))
                            % mark the prediction to be eliminated later
                            predictionMask(k) = true;
                            % eliminate index of prediction from table
                            outResults.resultsTable.prediction(indexes(k)) = 0;
                            continue
                        else
                            aggrLabel = aggrLabel.(S);
                        end
                    end
                    if ~predictionMask(k)
                        if isfield(pResults,'full')
                            oResults.full=...
                            aggregate(pResults.full,aggrLabel,aggregation{1});
                        end
                        % aggregate prediction labels corresponding to
                        % crossvalidation
                        if isfield(pResults, 'cv')
                            if length(pResults.cv)==1                                        
                                oResults.cv=...
                                    aggregate(pResults.cv,aggrLabel, aggregation{1});                                        
                            else
                                oResults.cv=cell(1,length(pResults.cv));
                                for i=1:length(pResults.cv)
                                    oResults.cv{i}=...
                                        aggregate(pResults.cv{i},aggrLabel, aggregation{1});
                                end
                            end
                        end
                        outResults.predictionResults{k}=oResults;
                    end
                end
            end
            % Remove predictions that could not be aggregated
            outResults.predictionResults = outResults.predictionResults(~predictionMask);
            % array containing indexes that remain with predictions
            newIndexes = find(~predictionMask);
            % asign new indexes to prediction column
            for i = 1:length(newIndexes)
                outResults.resultsTable.prediction(outResults.resultsTable.prediction == newIndexes(i)) = i;
            end

            % The same for aggregation of groundtruth
            
            %aggregate groundtruth labels            
            nDataPartition = length(obj.dataPartition);
            partitionMask = false(1,nDataPartition);
            outResults.dataPartition = cell(1,nDataPartition);
            propList = properties(MSDataPartition);
            for i=1:nDataPartition
                if ~isempty(obj.dataPartition{i})
                    if ischar(S)
                        if ~isempty(obj.dataPartition{i}.(S))
                            aggrLabel = obj.dataPartition{i}.(S);
                        else
                            partitionMask(i) = true;
                            outResults.resultsTable.data_partition(outResults.resultsTable.data_partition == i) = 0;
                            continue
                        end
                    end
                    outResults.dataPartition{i} = MSDataPartition;
                    for j=1:length(propList)
                        propValue = obj.dataPartition{i}.(propList{j});
                        if isa(propValue, 'MSLabelData')
                            outResults.dataPartition{i}.(propList{j}) =...
                            propValue.aggregate(aggrLabel,aggregation{2});
                        end
                    end
                end
            end
            % Remove data partitions that could not be aggregated
            outResults.dataPartition = outResults.dataPartition(~partitionMask);
            % array containing indexes that remain with predictions
            newIndexes = find(~partitionMask);
            % asign new indexes to prediction column
            for i = 1:length(newIndexes)
                outResults.resultsTable.data_partition(outResults.resultsTable.data_partition == newIndexes(i)) = i;
            end
        end
        
        function classifResultLabel = score2Label(obj,fun)
            % Converts the score results into label results
            % classifResultLabel=obj.score2Label % uses 'max' probability criterium
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
            
            if ~strcmp(obj.predictionType,'scores')
                disp('This method only applies for score prediction results')
                classifResultLabel=[];
                return
            end
            classifResultLabel=obj.copy;
            narginchk(1,2)
            if nargin==1
                classifResultLabel.transformPrediction(@score2Label);
            else
                try
                    classifResultLabel.transformPrediction(@score2Label,fun);
                catch
                    error(['The input parameter <fun> is not an appropriate'...
                          ' function to convert scores to label'])
                end
            end
            classifResultLabel.predictionType='label';
        end
        
        function aggregationROC=ROC4Aggregation(obj, indexes, S, cvIndex, once)
            % Returns an MSROCCurve object containing the ROC curves for a
            % specified aggregation.
            % aggregationROC=obj.ROC4Aggregation(preprocId,...
            %        featureExtId,mapModId,modelId, S, allData, plotFlag)
            % INPUT
            %   preprocId: string or index identifying the preprocessing
            %   featureExtId: string or index identifying the feature
            %       extraction method
            %   mapModId: string or index identifying the map modification
            %   modelId: string or index identifying the classification
            %       method
            %   S: MSLableData object indicating the level of aggregation
            %   allData: (optional) logical indicating whether the
            %       predictions correspond to the training and testing 
            %       using the whole data (default=false)
            %   plotFlag: (optional) logical indicating whether a plot with
            %       the performance measures is shown (default=true)
            % OUTPUT
            %   aggregationROC: MSROCCurve containing the ROC curves for a
            %       given aggregation
            
            if ~strcmp(obj.predictionType,'label')
                disp('This function is only available for label prediction')
                aggregationROC=[];
                return
            end
            %input validation
            narginchk(3,5)
            if nargin < 5
                once = false;
            elseif ischar(once) && strcmp(once,'once')
                once = true;
            else
                assertLogicalScalar(once, 'once');
            end
            if nargin<4 || isempty(cvIndex)
                cvIndex=0;
            else
                cvIndex=obj.validateCVIndex(cvIndex);
            end
            if ischar(S)
                if ~any(strcmp(S,{'individuals','samples','sampleSets'}))
                    error('Unknown aggregation type')
                end
            else
                aggrLabel = S;
            end
            
            % obtain predictions
            prediction=obj.prediction(indexes,cvIndex,once);
            if isempty(prediction)
                disp(['No predictions were stored for the input combinat'...
                    'ions. Therefore it is not possible to generate ROC curves'])
                aggregationROC = [];
                return
            end
            % obtain data partition
            [dPartition, indexPartition] = obj.findObject(indexes, 'data_partition');


            % obtain number of features vector
            [numberFeatures, indexFeatures] = obj.findObject(indexes, 'number_features');
            nPredictions=size(prediction,1);
                       
            aggregationROC=cell(nPredictions,1);

            for i=1:nPredictions
                if ~isempty(prediction{i})
                    % Compute ROC curves for each class prediction depending on
                    % aggregation threshold:
                    partitioni = dPartition{indexPartition(i)};
                    groundTruth=partitioni.classes;
                    nNumFeatures=length(numberFeatures{indexFeatures(i)});
                    if ischar(S)
                        if ~any(strcmp(S,fieldnames(partitioni)))...
                                ||isempty(partitioni.(S))
                            error(['The property <' S '> does not exist'])
                        end
                        aggrLabel = partitioni.(S);
                    end
                    % Convert data to single prediction if an external
                    % partition was use for the classification task and
                    % cvIndex==0
                    if iscell(prediction{i}) 
                        N = numel(prediction{i});
                        szPred = size(prediction{i}{1}.data);
                        pData = zeros(N*szPred(1),szPred(2));
                        gtData = zeros(N*szPred(1),1);
                        aggData = zeros(N*szPred(1),1);
                        for k=1:N
                            pData((szPred(1)*(k-1)+1):szPred(1)*k,:)=prediction{i}{k}.data;
                            gtData((szPred(1)*(k-1)+1):szPred(1)*k,:)=groundTruth.data;
                            aggData((szPred(1)*(k-1)+1):szPred(1)*k,:)=aggrLabel.data;
                        end
                        prediction{i}=MSLabelData(prediction{i}{1}.labels,pData);
                        gt = MSLabelData(groundTruth.labels,gtData);
                        sAgg = MSLabelData(aggrLabel.labels,aggData);
                    else
                        gt = groundTruth;
                        sAgg = aggrLabel;
                    end

                    aggregationClasses=gt.aggregate(sAgg,'unique');
                    if any(isnan(aggregationClasses.data))
                        error(['The aggregated groundtruth cannot be '...
                              'uniquely determined from the input aggregation label'])
                    end
                    % number of classes
                    nClasses=groundTruth.numLabels;

                    aggregationPrediction = zeros(sAgg.numLabels, prediction{i}.numLabels, ...
                                                              nNumFeatures);
                    for k = 1:nNumFeatures
                      % Item mask
                      M = prediction{i}.data(:,k) > 0;
                      % compute maximum from S
                      maxLabelUsed=max(sAgg.data(M,1));
                      % Contingency table
                      aggregationPrediction(1:maxLabelUsed,:,k) = ...
                                accumarray([sAgg.data(M,1), prediction{i}.data(M,k)], 1);
                    end

                    aggregationROC{i}=cell(1,nClasses);
                    % create an MSROCCurve for each class
                    M = aggregationClasses.data(:,1) > 0;
                    for posLabel=1:nClasses                
                        predScore = squeeze(aggregationPrediction(M, posLabel, :) ./ ...
                                            sum(aggregationPrediction(M,:,:), 2));
                        aggregationROC{i}{posLabel} = MSROCCurve(predScore, ...
                                              aggregationClasses.data(M,1) == posLabel,...
                                                    dPartition{indexPartition(i)}.classes.labels{posLabel});
                    end
    %                 %plot curves
    %                 figure(i)
    %                 for j=1:nClasses
    %                     subplot(1,nClasses,j)
    %                     aggregationROC{i}{j}.plotCurves(numberFeatures{indexFeatures(i)});                    
    %                 end
                end
            end
            if numel(aggregationROC)==1 && once
                aggregationROC=aggregationROC{1};
            end
        end       
        
        function setDataPartition (obj, dataPartition)
            if ~(isempty(obj.dataPartition) && isempty(obj.cvSegments))
                error(['Setting a dataPartition is only possible if '...
                       'no data partition or cv-segments are previously defined'])
            end
            % add data partition
            obj.dataPartition = {dataPartition};
            t = array2table(ones(obj.nCombinations,1));
            t.Properties.VariableNames = {'data_partition'};
            obj.resultsTable = horzcat(obj.resultsTable, t);
        end
        
        function names = selectorsFromIndexes (obj, indexes)
            % returns a list with strings identifying the sequence of
            % combinations for each index
            % INPUT
            %   indexes: vector of doubles representing the row of
            %   combinations in resultsTable
            % OUTPUT
            %   names: list of strings
           
            % check input
            if ~isnumeric(indexes) || ~all(indexes <= obj.nCombinations & indexes > 0)
                error(['The input <indexes> must be a numerical vector with indexes',...
                    ' in the range of the number of combinations'])
            end
            names = cell(1, length(indexes));
            t = obj.resultsTable(indexes,:);
            sNames=obj.selectorNames;
            names='';
            for j=1:length(sNames)
                names = strcat(names, '/', t.(sNames{j}));
            end
        end
    end
    
    methods (Access=protected)
        
        function validateDataPartition(~, dPartition)
            if ~isempty(dPartition) && ~isa(dPartition, 'MSDataPartition')
                error('Invalid input data partition')
            end
        end
        
        function indexes = validateIndexes (obj,indexes, once)
            narginchk(2,3)
            if nargin < 3
                once = false;
            elseif ischar(once) && strcmp(once,'once')
                once = true;
            else
                assertLogicalScalar(once, 'once');
            end
            if isempty(indexes)
                error('input <indexes> must not be empty')
            end 
            if iscellstr(indexes) 
                if length(indexes) > obj.nSelectors
                    error('Too many indexes as input')
                else
                    indexes=obj.select(indexes{:});
                end               
            elseif ~(isnumeric(indexes) && isvector(indexes) &&...
                    all(floor(indexes)==indexes) &&...
                    all(indexes>=1 & indexes <=obj.nCombinations))
                error(['Input <indexes> is incorrect. Make sure that '...
                    'the number of indexes is at most %d and that each '...
                    'index is a positive integer not greater than %d']...
                    ,obj.nSelectors, obj.nCombinations)
            end
            if once && ~isempty(indexes)
                indexes=indexes(1);
            end
        end         
        
        function addSelectorsFrom(obj1,obj2)
                        
            % get selector column names
            varNames2=obj2.resultsTable.Properties.VariableNames;
            varNames1=obj1.resultsTable(:,1:obj1.nSelectors).Properties.VariableNames;
            
            for i=1:obj2.nSelectors
                index=find(strcmp(varNames2{i},varNames1),1);
                if isempty(index) % the column does not exist in t1
                    obj1.addSelector(varNames2{i});
                end
            end
        end
        
        function [list, propertyName] = listFromColumn(obj,columnName)
            switch columnName                
                case 'preprocessing'
                    list = obj.preprocessings;
                    propertyName = 'preprocessings';
                case 'map'
                    list = obj.maps;
                    propertyName = 'maps';
                case 'map_modified'  
                    list = obj.mapsModified;
                    propertyName = 'mapsModified';
                case 'feature_data'  
                    list = obj.featureData;
                    propertyName = 'featureData';
                case 'model' 
                    list = obj.classifModels;
                    propertyName = 'classifModels';
                case 'data_partition'
                    list=obj.dataPartition;
                    propertyName = 'dataPartition';
                case 'cv_segment'
                    list=obj.cvSegments;
                    propertyName = 'cvSegments';
                case 'prediction'
                    list=obj.predictionResults;
                    propertyName = 'predictionResults';
                case 'number_features'
                    list=obj.numFeatures;
                    propertyName = 'numFeatures';
                otherwise
                    disp('Unknown type of object')
                    list = [];
            end
        end
        
        function addEmptyStorage(obj, columnName)
            c = zeros(obj.nCombinations,1);
            t = array2table(c);
            t.Properties.VariableNames = {columnName};
            obj.resultsTable=[obj.resultsTable t];
        end
        
        function transformPrediction(obj,fun,varargin)
            % Returns an MSClassificationResults object which modifies the
            % predictions contained in obj.predictionResults according to
            % the input function handle. The returned object has the same
            % property values as the input obj.
            % INPUT
            %   fun: function handle that receives an
            %   MSPredictionData and possibly other input parameters
            %   (specified in varargin) and returns a modified prediction
            %   varargin: The rest of input parameters of fun
            
            %aggregate prediction labels
            for k=1:length(obj.predictionResults)
                % if trained with all data aggregate the
                % prediction labels corresponding to that
                % training
                pResults=obj.predictionResults{k};                    
                if isfield(pResults,'full')
                    oResults.full=...
                    fun(pResults.full,varargin{:});
                end
                if isfield(pResults,'cv')
                    % aggregate prediction labels corresponding to
                    % crossvalidation
                    if length(pResults.cv)==1                                        
                        oResults.cv=...
                            fun(pResults.cv, varargin{:});                                        
                    else
                        oResults.cv=cell(1,length(pResults.cv));
                        for i=1:length(pResults.cv)
                            oResults.cv{i}=...
                                fun(pResults.cv{i}, varargin{:});
                        end
                    end
                end
                obj.predictionResults{k}=oResults;
            end
        end
        
    end
       
    methods   
      function objNew = copyUSETHIS ( obj )
        % Copy method leads to shallow copies for predictionResult fields
        % even though MSPredictionData < MSData < matlab.mixin.Copyable;
        % also there is no way to redefine the copy method here since the
        % copy method is sealed in matlab.mixin.Copyable
        % As a dirty workaround this copy method has been implemented to
        % ensure changes on the copied object won't effect the original
        % object
        objNew = obj.copy;
        for i = 1:objNew.nCombinations
          fieldList = fields( objNew.predictionResults{i} );
          for k = 1:numel(fieldList)
            try % No copy method for cv fields
              objNew.predictionResults{i}.(fieldList{k}) = ...
                objNew.predictionResults{i}.(fieldList{k}).copy;
            end
          end
        end        
      end
      
      function resultPart = reduce( obj, mask )
        % Reduce results to mask of items (logical vector), effecting the
        % predictionResults and the dataPartitions
        resultPart = obj.copyUSETHIS;
        for i = 1:resultPart.nCombinations
          fieldList = fields( resultPart.predictionResults{i} );
          for k = 1:numel(fieldList)
            resultPart.predictionResults{i}.(fieldList{k}).reduce( mask );
          end
        end
        for i = 1:numel(obj.dataPartition)
          resultPart.dataPartition{i} = ...
            resultPart.dataPartition{i}.reduce( mask );
        end
      end

      function resultParts = split( obj, labelData )
        % Cell Array of restrictions (.reduce) to multi-columned 
        % (!) labelData object
        resultParts = cell( labelData.dataLength, 1 );
        for i = 1:labelData.dataLength
          mask = logical( labelData.data(:,i)==2 );
          resultParts{i} = obj.reduce(mask);
        end
      end
    end
    
    methods (Static)
        
        function plotPerformance(performanceResults, groundTruth)
            %Plots performance results
            % INPUT
            %   performanceResults: struct with the following fields:
            %     -accuracy: list of positive doubles
            %     -kappa: list of positive doubles
            %     -sensitivity: matrix of positive doubles
            %     -specificity: matrix of positive doubles
            %     -avgSensitivity: list of positive doubles
            %     -avgSpecificity: list of positive doubles
            %     -numFeatures: list of number of features used for
            %     classification
            %   groundTruth: MSLabelData containing ground truth
            
            %Plot accuracy and kappa
            figure('name','General classification performance')%General performance
            plot(performanceResults.numFeatures,performanceResults.accuracy,'b',...
                 performanceResults.numFeatures,performanceResults.kappa,'r')
             title('General classification performance')
             legend('accuracy','kappa')
             nClasses=groundTruth.numLabels;
             figure('name','Performance measures per class');
             %Plot sensitivity and specificity per class
             for i=1:nClasses
                 className=groundTruth.labels{i};
                 subplot(1,nClasses,i)
                 plot(performanceResults.numFeatures,performanceResults.sensitivity(i,:))
                 hold on
                 plot(performanceResults.numFeatures,performanceResults.specificity(i,:))
                 title(className)
                 legend('sensitivity','specificity')
             end    
        end
        function performance=performanceFromLabel(prediction, groundTruth)
            % Return a performance struct given a prediction MSPredictionData
            % object.
            % INPUT
            %   predictionLabel: MSLabelData object containing the
            %       prediction for a classification task.
            %   performance: struct with the following fields:
            %      -accuracy: list of positive doubles
            %      -kappa: list of positive doubles
            %      -sensitivity: matrix of positive doubles
            %      -specificity: matrix of positive doubles
            %      -avgSensitivity: list of positive doubles
            %      -avgSpecificity: list of positive doubles
            %      groundTruth: MSLabelData containing groundTruth
            
            % obtain Confusion Matrix
            [~,confusionMatrix]=prediction.validate(groundTruth);
            % compute performance
            performance = MSConfusionMeasure(confusionMatrix);
            if isa(prediction, 'MSScoreData')
              try
                performance.AUC = prediction.ROC_AUC(groundTruth);
              catch
                performance.AUC = [];
              end
            end
            % plot results
            % MSClassificationResults.plotPerformance(performance, groundTruth);
        end
    end
    
    methods (Static)%,Access=protected)
        
        function mI = removeIndexes(indexes, inputMask)
            % eliminates from the list of indexes the positions which are
            % false in inputMask, shifting the indexes
            N = numel(inputMask);
            cumMask = cumsum(inputMask);
            mI=indexes;
            for i=1:N
                mask=indexes==i;
                mI(mask)=mI(mask)-i+cumMask(i);
            end
            % find indexes for which the mask is false (and should be removed)
            ind = find(~inputMask);
            for i=1:length(ind)
                mI(indexes==ind(i))=0;
            end
        end
        
        function cvIndex = validateCVIndex(cvIndex)
            % Validates the cvIndex for the search and converts it if
            % necessary to the corresponding integer: -1 for working with
            % the classification using the whole data, 0 for the
            % classification per cv and an integer if a specific cv is
            % required.
            if ischar(cvIndex)
                if strcmp(cvIndex,'full')
                    cvIndex=-1;
                elseif strcmp(cvIndex, 'cv')
                    cvIndex=0;
                else
                    error('Unknown <cvIndex> input type')
                end
            % includes the case where cvIndex is an array of positive
            % integers
            elseif (~(isscalar(cvIndex) && isnumeric(cvIndex) &&...
                    floor(cvIndex)==cvIndex && cvIndex>=-1))&&~ispositiveNaturalArray(cvIndex)
                error(['The input <cvIndex> must be a numerical '...
                  'scalar larger than -1'])
            end
        end
    end
end

function assertLogicalScalar(flag,varName)
    % Checks whether a variable, with name 'varName' in the
    % scope where the function is called is a logical scalar value. Throws an
    % error in case it doesn't
    if ~islogical(flag)||~isscalar(flag)
        error(['The Input <' varName '> must be a logical scalar'])
    end
end

function mergedPrediction = mergePredictions(oldPrediction,newPrediction)
    pFields = fields(oldPrediction);
    if ~all( strcmp(pFields,fields(newPrediction)) )
        error( 'Prediction fields do not match, possibly .full or .cv only existing once' );
    end
    mergedPrediction = oldPrediction;
    for i = 1:numel( pFields )
      oldLabelObj = oldPrediction.(pFields{i});
      newLabelObj = newPrediction.(pFields{i});
      mergedPrediction.(pFields{i}) = concatenatePredictionObjects( oldLabelObj, newLabelObj );
    end
end

function concatenatedPredictionObj = concatenatePredictionObjects( oldLabelObj, newLabelObj )
% Can concatenate MSLabelData as well as MSScoreData - position information
% will be deleted.
      labelsOld = oldLabelObj.labels;
      labelsNew = newLabelObj.labels;
      if ~all(strcmp(labelsOld,labelsNew))
        error('Labels do not match');
        % If only the order is different one could implement a resorting
        % This would need to take into account MSScoreData as well!
      end
      concatenatedPredictionObj = oldLabelObj.copy;
      concatenatedPredictionObj.data = vertcat( oldLabelObj.data, newLabelObj.data );
      concatenatedPredictionObj.setPositions([]);
end

function combinedDPartition = concatDPartitions( oldDPartition, newDPartition )
    if ~( strcmp(oldDPartition.specificLabels,'classes') && ...
        strcmp(newDPartition.specificLabels,'classes') )
      error('.specificLabels property in DPartition objects is not always pointing to "classes" ');
      % If this is desired this function needs to concatenate the respective 
      % label objects that are pointed to and the new DataPartition object
      % needs to decide for one label object ( e.g. classes ) that builds 
      % the concatenation
    end
    combinedDPartition = MSDataPartition();
    combinedDPartition.classes = concatenatePredictionObjects( oldDPartition.classes, newDPartition.classes );   
end


function bool = ispositiveNaturalArray(cvIndex)
bool = isvector(cvIndex) && isnumeric(cvIndex) &&...
                    all(floor(cvIndex)==cvIndex) && all(cvIndex>=1);
end
