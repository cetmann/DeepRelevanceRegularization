classdef MSClassificationTask < matlab.mixin.Copyable
 % This class represents different combinations of steps in the 
 % classification pipeline with different preprocessing sequences 
 % (normalization, mean filtering, spacial denoising, ...), feature 
 % extraction methods (Roc, PCA, NMF, ...), map modifications
 % (basis sorting, basis baseline removal, ...), and classification models 
 % (LDA, SVM, ...). Through the apply method the classification for all 
 % combinations is conducted and stored in an MSClassificationResults
 % object.
 %
 % Properties:
 %   preprocTable:  preprocessing table
 %   featureExtTable: feature extraction table
 %   mapModifTable: map modification table
 %   classifTable: classification model table
 %   numFeatures: number of features array
 %   processSteps: bool array identifying which process steps are to be
 %      conducted
 %   isRestricted: function with logical output indicating whether a sequence
 %      of processes should not be conducted (to avoid unnecessary
 %      combinations
 %
 % Methods
 %   apply: Conducts the classification task and return results as an
 %   MSClassificationResults object
 %   removeCombinations: restricts combinations that should not be
 %                         conducted
 %
 % MSClassificationTask uses the handle semantic, i.e. when 
 % assigning an object of this class to a variable, only a reference to
 % the original object is copied. Use the copy method to create a deep copy.
    properties (SetAccess=protected)
        preprocTable=[];
        featureExtTable=[];
        mapModifTable=[];
        classifTable=[];
        numFeatures=inf;
        processSteps=[true true true true];
    end    
    properties (Dependent)
        nPP;
        nFE;
        nMM;
        nC;
        listCombinations;
        listValidCombinations;
    end
    
    properties (SetAccess = protected)
        validCombinations
    end
    methods               
        function obj=MSClassificationTask(preprocTable, featureExtTable,...
                                               mapModifTable, classifTable, numFeatures)
            %   Tables (MSTable with properties <string>, containing a cell array of 
            %   strings which identify the processing step and <object>, which is the 
            %   object which applies the processing):
            %   preprocTable: <object> field contains an array of 
            %                   MSPreprocessingTrigger objects.
            %   featureExtTable: <object> field contains an MSFeatureExtraction
            %   mapModifTable: <object> field contains an array of MSMapModifier
            %                   objects.
            %   classifTable: <object> field contains an MSClassifier object

                                           
            narginchk(4,5)
            if nargin>4 && ~isempty(numFeatures)
                obj.numFeatures=numFeatures;
            end
            if ~isempty(classifTable)
                obj.classifTable=classifTable; 
            else
                error('At least a classification method must be defined')
            end
            if isempty(mapModifTable)
                mapModifTable=MSTable;
                mapModifTable.add('none', MSDummyMapModifier);
                obj.processSteps(3)=false;               
            end
            if isempty(featureExtTable)
                featureExtTable=MSTable;
                featureExtTable.add('none', MSDummyFeatureExtraction); 
                obj.processSteps(2)=false;
            end
            if isempty(preprocTable)
                preprocTable=MSTable;
                preprocTable.add('none', MSDummyTrigger);
                obj.processSteps(1)=false;
            end
            obj.preprocTable=preprocTable;
            obj.featureExtTable=featureExtTable;
            obj.mapModifTable=mapModifTable;
            
            obj.validCombinations = true(1,obj.nPP*obj.nFE*obj.nMM*obj.nC);
        end
        
        function restoreCombinations(obj)
            % Removes any restrictions over the executed combinations. After
            % the execution of this methods all the possible combinations 
            % (according to the tables used to create the classification
            % tasks) are included in the classification study
            obj.validCombinations = true(size(obj.validCombinations));
        end
        function removeCombinations (obj, varargin)
            % Removs the process combinations that should not be
            % considered according to the input patterns. The removed
            % combinations are those matching all of the given patterns
            % (logical AND applied to each filter)
            % obj.removeCombinations(input1, input2,...)
            % INPUT
            %   varargin: strings containing the pattern for the process ID
            %       e.g.,
            %       obj.removeCombinations('tic','','sort*')
            %       If the empty string ('') or the empty symbol ([]) are
            %       entered then all the id's for the corresponding
            %       cathegories are considered.            
                      
            processes = find(obj.processSteps);
            nProcesses = sum(obj.processSteps);
            narginchk(2,nProcesses+1);
            if ~iscellstr(varargin)
                error('All input arguments must be strings')
            end
            %init match indexes
            match = cell(1,4);
            for i = 1:4
                string = obj.getProcessString(i);
                match{i} = 1:numel(string);
            end
            
            % search restricted patterns
            for i = 1:length(varargin)
                string = obj.getProcessString(processes(i));
                if ~isempty(varargin{i})
                    restrictedIndex = (regexp(string,varargin{i},'once'));
                    restrictedIndex = find(~cellfun(@isempty,restrictedIndex));
                    if isempty(restrictedIndex)
                        return
                    end
                    match{processes(i)} = restrictedIndex;
                else
                    match{processes(i)} = 1:numel(string);
                end                
            end
            
            % asign false in the corresponding position for each
            % combination
            for i=match{1}
                s1 = (i-1)*obj.nFE*obj.nMM*obj.nC;
                for j=match{2}
                    s2 = s1 + (j-1)*obj.nMM*obj.nC;
                    for k=match{3}
                        s3 = s2 + (k-1)*obj.nC;
                        obj.validCombinations(s3 + match{4}) = false;
                    end
                end
            end
        end
        
        function [results]=apply(obj, data, dataPartition, kFold, varargin)
        % Performs a classification task for different combinations of
        % classification steps (preprocessing, feature extraction, map 
        % modifications and classification models)
        % INPUT
        %   data: MSMaldiData used for executing the classification task.
        %                 As a second option it can be a cell array with two string
        %                 elements, the first containing the name of the mat file
        %                 where the maldi data is stored and the second one
        %                 specifying the name of the variable
        %   dataPartition: MSDataPartion containing relevant labels for
        %                  the classification task

        %   kFold: Integer, matrix or MSLabelData indicating the
        %       kFold-crossvalidation
        %   varargin (Opcional input as pairs name-value)
        %   -trainAllData: bool indicating whether all data should be used for
        %                   training, independently from the crossvalidation.
        %   -featureExtrPerSegment: bool  whether feature extraction is
        %                   performed once with all the data or with each subset 
        %                   of training data for each crossvalidation segment
        %   
        %   -saveMap:   bool indicating whether generated maps are stored
        %   -saveFeatureData:  bool indicating whether generated feature data are stored
        %   -saveModels: bool indicating whether generated models are stored
        %   -predictionType: string indicating whether scores for fussy
        %        classification, labels for deterministic classification, or no
        %        predictions should be computed ('scores', 'label'[default], 'none')
        %  OUTPUT
        %   results: MSClassificationResults object containing properties and
        %             methods to evaluate the classification performance
        %
        %   In ClassifTest can be found an example where the input variables are
        %   initialized for a Lung tumor vs Pancreas tumor classification problem
        
        %----------------input validation----------------------------------      
        sourceDataOnDisk=false;
        if iscellstr(data)
            if ~(isvector(data)&&length(data)==2)
                error('The cell array data must contain only two strings')
            end
            try
                s=load(data{1},data{2});
                sourceData=s.(data{2});
            catch
                error(['Error by loading the maldi data. Make sure that the'...
                       'data path and the variable name are correct'])
            end  
            sourceDataOnDisk=true;
        elseif ~isa(data, 'MSMaldiData')
            error(['source data must be either an MSMaldiData object or a cell'...
                  'array of two strings'])
        else
            sourceData=data;
        end
        sourceData.assert;
        if ~isa(dataPartition, 'MSDataPartition')
            error('dataPartition must be an MSDataPartition object')
        end
        
        dataPartition.assert;  

        
        %------define other optional parameters as name-value pair-----------------
        params=inputParser;
        params.addParameter('trainAllData', true, @islogical);
        params.addParameter('featureExtrPerSegment', true, @islogical);
        params.addParameter('saveMaps', true, @islogical);
        params.addParameter('saveFeatureData', true, @islogical);
        params.addParameter('saveModels', true, @islogical);
        params.addParameter('predictionType', 'label', @obj.checkPredictionType);
        params.parse(varargin{:});
        %--------------------------------------------------------------------------
        
        autoPartition=false;
        if~isempty(kFold)
            if isa(kFold,'MSLabelData')
                MSClassificationTask.validatePartition(kFold.data, sourceData.numItems);            
            elseif ismatrix(kFold) && ~isscalar(kFold)
                MSClassificationTask.validatePartition(kFold, sourceData.numItems);
                kFold=MSLabelData(cellstr(num2str((1:size(kFold,2))')),kFold);
            elseif ~MSClassificationTask.isScalarGE(kFold,2)
                error('kFold must be an integer > 1')
            else
                kFold=floor(kFold);
                autoPartition=true;            
            end    
        elseif params.Results.featureExtrPerSegment 
%             choice = questdlg(['It is not possible to extract features per'...
%                 ' cv-segment, as no kFold number or partition was specified'], ...
%                 'Alert!', ...
%                 'Continue using the whole dataset','Abort study','Continue using the whole dataset');
%                 % Handle response
%             switch choice
%             case 'Abort study'
%                 disp('Study aborted by user')
%                 return
%             end
          warning(['It is not possible to extract features per'...
                 ' cv-segment, as no kFold number or partition was specified', ...
                 'Continuing using the whole dataset']);
        end
        % for the ease of notation new variables are defined
        npp=obj.nPP;
        nfe=obj.nFE;
        nmm=obj.nMM;
        nc=obj.nC;
        nTotalM=npp*nfe;
        nTotalMM=nTotalM*nmm;
        nTotalC=nTotalMM*nc;
        featurePerSegment=params.Results.featureExtrPerSegment && ~isempty(kFold);
        featureAllData=~featurePerSegment || params.Results.trainAllData;
        saveMaps=params.Results.saveMaps;
        saveFD=params.Results.saveFeatureData;
        saveModels=params.Results.saveModels;
        numClasses=dataPartition.classes.numLabels;
        computePrediction = (~strcmp(params.Results.predictionType,'none') || ...
                    ~isempty(kFold)) && ~(isempty(kFold) && ~params.Results.trainAllData);  
        %initialize output variable
        results=MSClassificationResults;
        results.predictionResults=cell(1,npp*nfe*nmm*nc);
        results.numFeatures={obj.numFeatures};
        results.resultsTable=table;        
        results.nSelectors=sum(obj.processSteps);
        stored=[obj.processSteps(1), saveMaps, saveMaps saveFD, ...
            saveModels, computePrediction, true, true, true];

        %%%%%%%%%%%%%%% defining cvSegment labels and related structures %%%%%%%%%%
        
        % cvSegments contains the partition of the data for the crossvalidation
        % currentTrainLabels and currentTestLabels are a cell array of labels for
        % training and testing for each crossvalidation partition.
        if ~isempty(kFold)
            [cvSegments, currentTrainLabels, currentTestLabels] = MSClassificationTask.MSPartition(dataPartition, kFold);
            numCVIterations = length(currentTrainLabels);
            %Store the cross validation partition
            results.cvSegments={cvSegments};
        else
            
        end
        %Store preprocessings
        if obj.processSteps(1)
            results.preprocessings=obj.preprocTable.object;
        end

        %Store data partition for evaluating performance and for
        %visualization
        results.dataPartition={dataPartition};
        
        if computePrediction
            results.predictionType=params.Results.predictionType;
        else
            results.predictionType='none';
        end
        if computePrediction && strcmp(results.predictionType,'none')
            disp(['The default label predictions will be computed despite the input '...
                ' prediction type because a cv-partition was specified'])
            results.predictionType='label';
        end

        %initialize maps, feature data and classif model arrays.

        if saveFD
            results.featureData=cell(1, nTotalMM);
        end
        if saveMaps
            results.maps=cell(1,nTotalM);
            results.mapsModified=cell(1,nTotalMM);
        end
        if saveModels 
            results.classifModels=cell(1,nTotalC);
        end
        % init storage indexes
        mapIndex=1;
        modifiedMapIndex=1;
        modelIndex=1;
        predictionIndex=1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %--------------pipeline implementation------------------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% preprocessing level (iPP)%%%%%%%%%%%%%%%%%%%%%
        firstPP = true;
        for iPP=1:npp 
            if ~obj.executeIteration(iPP)
                continue
            end
            disp(['Start preprocessing "' obj.preprocTable.string{iPP} '"...'])
            % prepare data object to work on
            if sourceDataOnDisk
                % source data is loaded from disk
                if firstPP
                    % data has been loaded already, use it without copying
                    trainingData = sourceData;
                    % remove reference of sourceData so that data can later
                    % be removed
                    clearvars sourceData;
                else
                    % re-load from disk
                    s=load(data{1},data{2});
                    trainingData=s.(data{2});
                end
            else
                % source data is in memory
                if obj.processSteps(1)
                    % work on a copy
                    trainingData=sourceData.copy;
                else
                    % no preprocessing specified, no need to copy
                    trainingData=sourceData;
                end
            end
            firstPP = false;
            %--perform current sequence of preprocessing
            currentPrepList=obj.preprocTable.object{iPP};
            try
            for i=1:length(currentPrepList)                
                currentPrepList(i).apply(trainingData);
            end
            catch e
                    obj.handleError(results, e, obj.preprocTable.string{iPP});
                    continue;
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% feature extraction level (iFE) %%%%%%%%%%%%%%%
            for iFE=1:nfe    
                if ~obj.executeIteration([iPP, iFE])
                    continue;
                end
                disp(['Start feature extraction "' obj.featureExtTable.string{iFE} '"...'])
                try
                projection=obj.featureExtTable.object{iFE};         
                catch e
                    errorStep = [obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE}];
                    obj.handleError(results, e, errorStep);
                    continue;
                end
                %--------------- Generate feature map and feature data-------------
                if featureAllData
                    try
                    % compute maps and feature data using all the data
                    if projection.supportsDataPartition
                        map = projection.createMap(trainingData, dataPartition);
                    else
                        map = projection.createMap(trainingData, dataPartition.classes);
                    end
                    if isa(map,'MSBasisMap')
                        decObj = map.decompositionObj;
                    else
                        decObj = [];
                    end
                    catch e
                        errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE}...
                                   ' during map generation'];
                        obj.handleError(results, e, errorStep);
                        continue;
                    end
                    %store resulting maps and features
                    if saveMaps
                        results.maps{mapIndex}.full=map;
                    end
                end
                if featurePerSegment
                    mapCV = cell(1, numCVIterations);
                    decObjCV = cell(1, numCVIterations);
                    % compute maps and feature data for CV-segment
                    for iCV = 1:numCVIterations 
                        try
                            if projection.supportsDataPartition
                                mapCV{iCV} = projection.createMap(trainingData, currentTrainLabels{iCV});
                            else
                                mapCV{iCV} = projection.createMap(trainingData, currentTrainLabels{iCV}.classes.data);
                            end
                            if isa(mapCV{iCV},'MSBasisMap')
                                decObjCV{iCV} = mapCV{iCV}.decompositionObj;
                            else
                                decObjCV{iCV} = [];
                            end
                        catch e
                            errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                       ', segment ' num2str(iCV) 'during map generation'];
                            obj.handleError(results, e, errorStep);
                            continue;
                        end
                    end
                    %store resulting maps and features
                        if saveMaps                
                            results.maps{mapIndex}.cv=mapCV;
                        end
                end
                mapIndex=mapIndex+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%% map modifier level (iMM) %%%%%%%%%%%%%%%%%%%%%%%
                for iMM=1:nmm
                    if ~obj.executeIteration([iPP, iFE, iMM])
                        continue;
                    end
                    disp(['Start map modification "' obj.mapModifTable.string{iMM} '"...'])
                    currentMapModList=obj.mapModifTable.object{iMM};
                    modified=true;
                    if featureAllData
                        try
                        % compute the sequence of map modifications
                        [mapMod, modified]=obj.applyMapModList(currentMapModList, ...
                                        map, trainingData, dataPartition.classes);
                        catch e
                            errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                       '+' obj.mapModifTable.string{iMM}];
                            obj.handleError(results, e, errorStep);
                            continue;
                        end
                        %it is possible that no modification was performed because it was not applicable
                        if modified 
                            try
                            featureDataMod=mapMod.map(trainingData);
                            catch e
                                disp(e.message)
                            end
                            %store resulting maps and features in arrays
                            if saveMaps
                                results.mapsModified{modifiedMapIndex}.full=mapMod;
                            end
                            if saveFD
                                results.featureData{modifiedMapIndex}.full=featureDataMod;
                            end
                        end
                    end
                    if featurePerSegment && modified
                        mapModCV=cell(1,numCVIterations);
                        featureDataModCV=cell(1,numCVIterations);
                        for iNM = 1:numCVIterations
                            try
                            [mapModCV{iNM}, modified]=obj.applyMapModList(currentMapModList,...
                                mapCV{iNM}, trainingData, currentTrainLabels{iNM});

                            % Generate feature data using the nmfMap
                            featureDataModCV{iNM} = mapModCV{iNM}.map(trainingData);
                            catch e
                                errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                           '+' obj.mapModifTable.string{iMM} num2str(iCV) ' during segment '...
                                           num2str(iNM)];
                                obj.handleError(results, e, errorStep);
                                continue;
                            end
                        end
                        %store maps and features per CV-segment
                        if saveMaps
                            results.mapsModified{modifiedMapIndex}.cv=mapModCV;
                        end
                        if saveFD
                            results.featureData{modifiedMapIndex}.cv=featureDataModCV;
                        end
                       
                    end
                    if modified 
                        modifiedMapIndex=modifiedMapIndex+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% classification level (iC) %%%%%%%%%%%%%%%%%%%
                        for iC = 1:nc 
                            if ~obj.executeIteration([iPP, iFE, iMM, iC]);
                                continue;
                            end
                            
                            classifier=obj.classifTable.object{iC};
                            
                            % check if decObjs are supported
                            % note: whole classifier row is skipped if any
                            % decObj is not supported, in general support
                            % should be same for all decObjs
                            if isa(classifier, 'MSDecompositionBasedClassifier')
                                if (params.Results.trainAllData || (~isempty(kFold) && ~featurePerSegment)) ...
                                    && ~classifier.supportsDecomposition(decObj)
                                    continue;
                                end
                                if ~isempty(kFold) && featurePerSegment
                                    dec_supp = true;
                                    for iCV=1:numCVIterations
                                        if ~classifier.supportsDecomposition(decObjCV{iCV})
                                            dec_supp = false;
                                            break;
                                        end
                                    end
                                    if ~dec_supp
                                        continue;
                                    end
                                end
                            end
                            
                            disp(['Start classification "' obj.classifTable.string{iC} '"...'])
                            
                            % scores required => classification supports
                            % scores
                            scoresSupp=classifier.supportScores || ~strcmp(results.predictionType,'scores');
                            computePrediction = stored(end-3) && scoresSupp;

                            % Add row to table
                            results.resultsTable=[results.resultsTable; table(obj.preprocTable.string(iPP),...
                              obj.featureExtTable.string(iFE), obj.mapModifTable.string(iMM),...
                              obj.classifTable.string(iC), iPP, mapIndex-1, ...
                              modifiedMapIndex-1, modifiedMapIndex-1,...
                              modelIndex, predictionIndex*scoresSupp,1,1,1)];
                          
                                                             % increment model index
                             modelIndex=modelIndex+1; 
                          
                            % Initialization of prediction structure if
                            % autoPartition flag removed (!)
                            % Prediction results are now always (!) saved
                            % in separate columns for each CVFold
                            if computePrediction && ~isempty(kFold)
                                %initialize the crossvalidation prediction labels
                                cvPrediction = obj.initPredictionData(currentTrainLabels,...
                                                        trainingData, length(obj.numFeatures),...
                                                        numClasses, results.predictionType,numCVIterations);
                                %Initialize cvSelfPrediction to
                                %to determine thresholds on the
                                %training folds later
                                cvSelfPrediction = obj.initPredictionData(currentTrainLabels,...
                                                        trainingData, length(obj.numFeatures),...
                                                        numClasses, results.predictionType, numCVIterations);
                            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of features level (k) %%%%%%%%%%%%%%%%%
                             for k=1:length(obj.numFeatures)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cross validation level (iCV) %%%%%%%%%%%%%%
                                if ~isempty(kFold)
                                  for iCV= 1:numCVIterations
                                    fprintf('Starting #%d out of %d CV iterations\n',iCV,numCVIterations);
                                    if featurePerSegment
                                        fd=featureDataModCV{iCV}.reduceFeatures(obj.numFeatures(k));
                                        do=decObjCV{iCV};
                                    else
                                        fd=featureDataMod.reduceFeatures(obj.numFeatures(k));
                                        do=decObj;
                                    end
                                    if isa(classifier, 'MSDecompositionBasedClassifier')
                                        classifier.setDecomposition(do);
                                    end
                                    try
                                    model = classifier.trainModel(fd,...
                                        currentTrainLabels{iCV}.classes);
                                    catch e
                                        errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                             '+' obj.mapModifTable.string{iMM} num2str(iCV)...
                                             '+' obj.classifTable.string{iC} ' during training with '...
                                             num2str(k) ' features at segment ' num2str(iCV)];
                                        obj.handleError(results, e, errorStep);
                                        continue;
                                    end
                                    if saveModels
                                        results.classifModels{modelIndex-1}.cv{k,iCV}=model;
                                    end
                                    if computePrediction
                                        % Apply model for predictions and scores if required
                                        testLabelMask = currentTestLabels{iCV}.classes.data > 0;
                                        trainLabelMask = currentTrainLabels{iCV}.classes.data > 0;
                                        % Select prediction object to be updated in this iteration
                                        % (Getting a handle)
                                        cvPredictionCurrent = cvPrediction{iCV};
                                        % Part self
                                        cvPredictionSelfCurrent = cvSelfPrediction{iCV};
                                       try
                                       if strcmp(results.predictionType,'scores')
                                           [~, scores] = model.classify(fd, testLabelMask); 
                                           cvPredictionCurrent.data(testLabelMask,numClasses*(k-1)...
                                             +1:numClasses*k) = scores.data(testLabelMask,:);
                                           % Part self
                                           [~, scoresSelf] = model.classify(fd, trainLabelMask);
                                           cvPredictionSelfCurrent.data(trainLabelMask,numClasses*(k-1)...
                                             +1:numClasses*k) = scoresSelf.data(trainLabelMask,:);
                                       else
                                           prediction = model.classify(fd, testLabelMask); 
                                           cvPredictionCurrent.data(testLabelMask,k) =...
                                                    prediction.data(testLabelMask,1);
                                           % Part self
                                           predictionSelf = model.classify(fd, trainLabelMask); 
                                           cvPredictionSelfCurrent.data(trainLabelMask,k) =...
                                                    predictionSelf.data(trainLabelMask,1);
                                       end
                                       catch e
                                            errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                                 '+' obj.mapModifTable.string{iMM} num2str(iCV)...
                                                 '+' obj.classifTable.string{iC} ' during predicting with '...
                                                 num2str(k) ' features at segment ' num2str(iCV)];
                                            obj.handleError(results, e, errorStep);
                                            continue;
                                       end
                                    end
                                  end
                                end
                             end

                             %Train and test model using all data if required
                             if params.Results.trainAllData
                                 if computePrediction
                                     if strcmp(results.predictionType,'scores')                             
                                         predictionAllData = ...
                                             MSScoreData(dataPartition.classes.labels,...
                                                    zeros(dataPartition.classes.numItems,...
                                             length(obj.numFeatures)*numClasses), trainingData);
                                     else
                                         predictionAllData = ...
                                            MSLabelData(dataPartition.classes.labels,...
                                                     zeros(dataPartition.classes.numItems,...
                                                    length(obj.numFeatures)), trainingData);
                                     end
                                 end
                                 if isa(classifier, 'MSDecompositionBasedClassifier')
                                     classifier.setDecomposition(decObj);
                                 end
                                 for k=1:length(obj.numFeatures)
                                     try
                                         fd = featureDataMod.reduceFeatures(obj.numFeatures(k));
                                         model = classifier.trainModel(fd,...
                                                dataPartition.classes);
                                     catch e
                                         errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                             '+' obj.mapModifTable.string{iMM}...
                                             '+' obj.classifTable.string{iC} ' during training with '...
                                             num2str(k) ' features'];
                                         obj.handleError(results, e, errorStep);
                                         continue
                                     end
                                     if saveModels
                                         results.classifModels{modelIndex-1}.full{k}=model;
                                     end
                                     if computePrediction
                                         try
                                         if strcmp(results.predictionType,'scores')
                                             [~, scores] =  model.classify(fd,...
                                                               dataPartition.classes.data>0); 
                                             predictionAllData.data(:,numClasses*(k-1)...
                                                                     +1:numClasses*k) = scores.data;
                                         else
                                             prediction = model.classify(fd,...
                                                        dataPartition.classes.data>0);
                                             predictionAllData.data(:,k)=prediction.data;
                                         end
                                         catch e
                                             errorStep=[obj.preprocTable.string{iPP} '+' obj.featureExtTable.string{iFE} ...
                                                 '+' obj.mapModifTable.string{iMM}...
                                                 '+' obj.classifTable.string{iC} ' during predicting with '...
                                                 num2str(k) ' features'];
                                             obj.handleError(results, e, errorStep);
                                             continue
                                         end
                                     end                                         
                                 end
                              end

                             %Build the structure "result" representing results for
                             %current parameter combination
                             if computePrediction
                                 if params.Results.trainAllData
                                     result.full=predictionAllData;
                                 end
                                 if ~isempty(kFold)
                                     if length(cvPrediction)==1
                                         cvPrediction = cvPrediction{1};
                                     end
                                     if length(cvSelfPrediction)==1
                                         cvSelfPrediction = cvSelfPrediction{1};
                                     end                               
                                     result.cv = cvPrediction;
                                     result.cvSelf = cvSelfPrediction;
                                 end
                                 %Store results
                                 results.predictionResults{predictionIndex}=result;                       
                                 predictionIndex=predictionIndex+1;                                     
                             end
                             disp(['Results for combination ' obj.printCombination([iPP,iFE,iMM,iC]) ' stored'])
                        end
                    end
                end
            end
        end
        if ~isempty(results.errorReport)
            results.errorReport=sprintf(results.errorReport);
        end
        % eliminate columns from the table that are not used and name
        % columns
        selector=[obj.processSteps stored];
        if isempty(results.resultsTable)
            warning('No results could be generated for the given combinations')
            results = [];
            return
        end
        results.resultsTable=results.resultsTable(:,selector);
        colNames={'PP', 'FE', 'MM', 'CM', 'preprocessing', 'map',...
                  'map_modified','feature_data', 'model',...
                  'prediction','data_partition','cv_segment', 'number_features'};
        results.resultsTable.Properties.VariableNames=colNames(selector);
        results.mainProcessSelectors = colNames(1:4);
        results.mainStorageIndexes = colNames(5:end);
        
        % eliminate empty endings of cell arrays containing objects
        if saveMaps
            results.maps = results.maps(1:mapIndex-1);
            results.mapsModified = results.mapsModified(1:modifiedMapIndex-1);
        end
        if saveFD
            results.featureData = results.featureData(1:modifiedMapIndex-1);
        end
        if saveModels
            results.classifModels = results.classifModels(1:modelIndex-1);
        end
        if predictionIndex > 1
            results.predictionResults = results.predictionResults(1:predictionIndex-1);
        else
            results.predictionResults = [];
        end
        disp('The classification study has been succesfully completed!!')
        end
    end
    
    methods %Property set and get
        
        %----------------<property set>------------------------------------
        function set.preprocTable(obj, ppT)
            MSClassificationTask.assertTable(ppT, 'MSPreprocessingTrigger', true);
            obj.preprocTable=ppT;
        end
        function set.featureExtTable(obj, feT)
            MSClassificationTask.assertTable(feT,'MSFeatureExtraction', false)
            obj.featureExtTable=feT;
        end
        function set.mapModifTable(obj, mmT)
            MSClassificationTask.assertTable(mmT, 'MSMapModifier', true);
            obj.mapModifTable=mmT;
        end
        function set.classifTable(obj, cT)
            MSClassificationTask.assertTable(cT, 'MSClassifier', false);
            obj.classifTable=cT;
        end
        function set.numFeatures(obj,numFeatures)
            if ~MSClassificationTask.isNumericVectorGE1(numFeatures)
                error('numFeatures must be a numeric positive vector')
            end
            obj.numFeatures=floor(numFeatures);
        end
        %---------------<\property set>------------------------------------
        
        %---------------<property get>-------------------------------------
        function value=get.nPP(obj)
            value=length(obj.preprocTable.string);            
        end
        function value=get.nFE(obj)
            %returns number of different feature extraction methods tested
            value=length(obj.featureExtTable.string);
        end
        function value=get.nMM(obj)
            %returns number of different map modifications tested
            value=length(obj.mapModifTable.string);
        end
        function value=get.nC(obj)
            %returns number of different classification models tested
            if isempty(obj.classifTable)
                value=[];
            else
                value=length(obj.classifTable.string);
            end
        end
        function value=get.listCombinations(obj)
          [a,b,c,d] = ndgrid( (1:obj.nC)', (1:obj.nMM)', ...
            (1:obj.nFE)', (1:obj.nPP)' );
          indexList = [d(:) c(:) b(:) a(:)];
          lIdx = size( indexList, 1 );
          value = cell(lIdx,1);
          for i = 1:lIdx
            value{i} = obj.printCombination(indexList(i,:));
          end 
        end
        function value=get.listValidCombinations(obj)
          value = obj.listCombinations;
          value = value(obj.validCombinations);
        end
        %------------------------------------------------------------------
        
        function string = printCombination(obj, indexes)
            string='';
            processTables={obj.preprocTable.string,obj.featureExtTable.string,...
                           obj.mapModifTable.string, obj.classifTable.string};
            for i=1:4
                if obj.processSteps(i)
                    string=strcat(string, processTables{i}{indexes(i)},'-');
                end
            end
            string=string(1:end-1);
        end
        
    end
    
    methods(Static, Access = private)
        
        function lbls = deselectLabels(lbls,mask)
            lbls = lbls.copy;
            lbls.data(mask) = 0;
        end
        
        function cvPrediction = initPredictionData(currentTrainLabels,...
                                                trainingData, nFeatures, numClasses, predictionType, nObj)
            cvPrediction=cell(1,nObj);  
            if strcmp(predictionType,'scores')
                ncol=nFeatures*numClasses;
                f=@MSScoreData;
            else
                ncol=nFeatures;
                f=@MSLabelData;
            end  
            for i=1:nObj
            cvPrediction{i} = ...
                f(currentTrainLabels{1}.classes.labels,...
                  zeros(currentTrainLabels{1}.numItems,ncol,'like',trainingData.data),trainingData);
            end           
        end
        
        function handleError(results, e, step)
            % Handle exception occurred while processing a method variant
            results.nErrors=results.nErrors+1;
            location = [strrep(e.stack(1).file, '\', '/') '#' num2str(e.stack(1).line)];
            msg=[e.message ' (' location ')'];
            errorMessage=['Error ' num2str(results.nErrors) ' occurred at '...
                          step '-: ' msg];
            results.errorReport=[results.errorReport '\n' errorMessage];
            % Temporarily disable warning backtrace
            pushState=warning('off', 'backtrace');
            warning(errorMessage);
            warning(pushState.state, pushState.identifier);
        end
    end
    
    methods(Static)
        
    function result = isScalarGE (X,num)
        result = isnumeric(X) && isscalar(X) && X >= num;
    end

    function result = isNumericVectorGE1 (X)
        result = isnumeric(X) && isvector(X) && ~any(X<1);
    end      

    function assertTable(X,className,allowEmptyObjects)
        % checks whether the given table X is a valid table with objects of class
        % className
        if ~isstruct(X)&&~isa(X,'MSTable')
            error('The table must be either a struct or an MSTable object')
        end
        fieldsX=fields(X);
        if ~(length(fieldsX)==2 && any(strcmp(fieldsX, 'string')) && any(strcmp(fieldsX, 'object')))
            error('The table struct must have exclusively the fields <string> and <object>')
        end
        if ~iscellstr(X.string)||~isvector(X.string)
            error('The field <string> must be a cell array of strings')
        end
        if ~iscell(X.object)||~isvector(X.object)
            error('The field <object> must be a cell array')
        end
        if length(X.string)~=length(X.object)
            error('Fields <string> and <object> must have the same length')
        end
        for i=1:length(X.object)
            if ~(isa(X.object{i}, className) || ...
                 (allowEmptyObjects && isempty(X.object{i})))
                error(['The field <object> must be a cell array of objects of type <' className '>'])
            end
        end
    end

    function [cvSegments, currentTrainLabels, currentTestLabels] = MSPartition(dataPartition, kFold) 
    % Given a MSDataPartition object dataPartition and an integer kFold
    % indicating number of crossvalidation segments returns:
    % -an MSLabel data cvSegments, indicating the data partition in cv-segments
    % -a cell array of MSLabelData currentTrainLabels where the kth column of the
    %   data indicates the training labels for the kth crossvalidation segment
    % -a cell array of MSLabelData currentTestLabels where the kth column of the
    %   data indicates the test labels for the kth crossvalidation segment

    %--define cross validation segments
    if isnumeric(kFold)
        if isempty(dataPartition.cvEntities)
            cvSegments = dataPartition.classes.kFold(kFold);
        else
            cvSegments = dataPartition.classes.kFold(kFold,dataPartition.cvEntities);
        end

        % Generate cross validation labels
        % For each of j=1..K cross validation segments, test set is selected by
        % cvSegments label == j, training set is the complement.
        numCVIterations = cvSegments.numLabels;
        testLabels = dataPartition.applyLabelFunction(@(lbls) lbls.spread(cvSegments));
        trainLabels = dataPartition.applyLabelFunction(@(lbls) testLabels.classes.complement(lbls));

        % Initialize containers for iteration results
        % cvResults = struct();
        currentTrainLabels=cell(numCVIterations, 1);
        currentTestLabels=cell(numCVIterations, 1);
        for iCV = 1:numCVIterations 
          % Select train and test labels for this iteration
          currentTrainLabels{iCV} = trainLabels.applyLabelFunction(@(lbls) lbls.select(iCV));
          currentTestLabels{iCV} = testLabels.applyLabelFunction(@(lbls) lbls.select(iCV));
        end
    else
        % If the partition was externally defined create train and test labels as the partition indicates
        [currentTrainLabels, currentTestLabels]=MSClassificationTask.splitFromPartition(dataPartition, kFold);
        cvSegments=kFold;
    end
    end
    
    function validatePartition(data, rows)
        if size(data,1)~=rows
            error('The number of data rows must coincide with the numItems of trainingData')
        end
        % verify that data contains only 0s, 1s and 2s
        for i=1:size(data,2)
            c=unique(data(:,i));
            if (length(c)==3 && c(1)~=0)||(~all(c(end-1:end)==[1; 2]))
                error(['Partition data columns must contain as non-zero entries '...
                    '1s and 2s for training and testing, respectively'])
            end      
        end
    end
    
    function [trainLabels, testLabels]=splitFromPartition(dataPartition, partition)
        nCV = partition.dataLength;
        trainLabels = cell(nCV,1);
        testLabels = cell(nCV,1);
        for iCV=1:nCV
            mask=(partition.data(:,iCV)==1);
            trainLabels{iCV}=dataPartition.applyLabelFunction(@(lbls) MSClassificationTask.deselectLabels(lbls,~mask));
            mask=(partition.data(:,iCV)==2);
            testLabels{iCV}=dataPartition.applyLabelFunction(@(lbls) MSClassificationTask.deselectLabels(lbls,~mask));
        end
    end
    
    function checkPredictionType(prediction)
        if ~(ischar(prediction)||any(strcmp(prediction,{'label','scores','none'})))
            error('Invalid prediction type. Valid types are: label, scores and none')
        end
    end

    function [modMap, modified]=applyMapModList(mapModList, map, data, labels)
    % applies the map modifications indicated in the list of MSMapModifier 
    % objects mapModList to the MSFeatureMap map.
    % The MSData data and MSLabels labels are used for certain map
    % modifications
    % if number of outputs is 2 returns a bool indicating whether it was possible
    % to perform at least one of the modifications

    modMap=map;
    if isa(modMap,'MSFeatureMapSequence')
        modMapPart2Change = modMap.listFM{end};
    else
        modMapPart2Change = modMap;
    end
    
    
    for i=1:length(mapModList)
        if mapModList(i).applicable(modMapPart2Change)
            modMapPart2Change = mapModList(i).apply(modMapPart2Change,data,labels);
        else
            if i>1
                warning(['One of the map modification combinations contains'...
                          ' modifications that can be applied to a feature'...
                          ' extraction method and others that cannot. The '...
                          ' corresponding results for this combination will'...
                          ' not be stored']);
            end
            modified=false;
            return;
        end
    end
    
    if isa(modMap,'MSFeatureMapSequence')
        modMap.listFM{end} = modMapPart2Change;
    else
        modMap = modMapPart2Change;
    end
    
    
    modified=true;
    end
    end
    
    methods(Access = protected)       
                
        function bool = executeIteration(obj,indexes)            
            nProcesses = [obj.nPP obj.nFE obj.nMM obj.nC];
            nIndex = length(indexes);
            d = prod(nProcesses(nIndex+1:end));
            mult = d;
            s = 0;
            for i= nIndex:-1:1
                s = s + (indexes(i)-1)*mult;
                mult = mult*nProcesses(i);
            end            
            bool = any(obj.validCombinations(s+1:s+d));
        end
        
        function string = getProcessString(obj, process, id)
            if nargin < 3
                id = 0;
            end
            switch process
                case 1
                    string = obj.preprocTable.string;
                case 2
                    string = obj.featureExtTable.string;
                case 3
                    string = obj.mapModifTable.string;
                case 4
                    string = obj.classifTable.string;
            end
            if id
                string = string{id};
            end
        end
    end
    
end
