classdef MSParameterSearch < matlab.mixin.Copyable
    % Parameter Search.
    %
    % This class implements a nested cross-validation for hyperparameter
    % optimization in order to compare different classification pipelines
    % to each other, where each pipeline has individually optimized
    % hyperparameters. A dataset is first partitioned into N different
    % subsets CV_1,...CV_N. Our goal is to make an aggregated prediction
    % for the whole dataset, where the prediction for each CV_i is based on
    % a model trained on Tr_i = union{CV_1,...,CV_(i-1),CV_(i+1),...,CV_N}.
    % A pre-specified set of parameters is then optimized via
    % cross-validation on each Tr_i.
    %
    % methods:
    %   MSParameterSearch: constructor. Takes as input trainingData
    %       (MSMaldiData), trainingPartition (MSDatapartition), parameters
    %       (cell of arrays of parameters), classificationPipeline
    %       (function handle of classification pipeline, which takes as
    %       input data,trainingLabels,testLabels,parameters,i,j).
    %   parameterOptimization: Searches the optimal parameters.
    %   conductStudy: If the parameters haven't been found yet, call
    %       parameterOptimization. Otherwise, apply the optimal parameters
    %       to the Tr_i.
    %   evaluate: Calculates the aggregated score.
    properties
        trainingData % The whole dataset.
        trainingPartition % A partition containing at least a field classes
                          % and a field cvEntities.
        parameters % A cell of parameter arrays.
        CVSubsegmentName % Field for the stratified parameter search.
                         % Defaults to 'individuals'.
        CVSubsegmentNumber % Number of inner CV segments.
        classificationPipeline % Function handle for the classification 
                               % pipeline. Should output a prediction.
        evaluationMetric % Defaults to 'avgSensitivity'.
        intermediateResultFile % If set, saves after every training.
        parameterSearchFinished % true or false.
        classNames % The class names.
        CV % Struct containing information about the CV.
        prediction % MSLabelData for the prediction.
        score % The predicted score according to the evaluationMetric.
        parameterChoiceModifier % Add this to the parameter index (within
                                % range) after finding the best parameter
                                % on the validation set, e.g. to choose
                                % a more robust model for the test set.
                                     
        trainLabels
        testLabels
    end
    properties (Access=private)
        kFoldMatrix
    end
    %% Plan for the parameter search
    
    methods
        function obj = MSParameterSearch(... 
                trainingData,...
                trainingPartition,...
                parameters,...
                classificationPipeline,...
                varargin) % Constructor
            
            obj.trainingData = trainingData;
            obj.trainingPartition = trainingPartition;
            obj.parameters = parameters;


            
            parser = inputParser;
            
            
            % @isstr does not work for input parsers.
            checkString = @(x) validateattributes(x,{'char'},{'nonempty'});
            
            addParameter(parser,'CVSubsegmentName','individuals',checkString)
            addParameter(parser,'CVSubsegmentNumber',2,@isnumeric);
            addParameter(parser,'evaluationMetric','avgSensitivity',checkString)
            addParameter(parser,'intermediateResultFile','',checkString)
            addParameter(parser,'parameterChoiceModifier',0,@isnumeric);
            
            parse(parser,varargin{:});
            
            obj.CVSubsegmentName = parser.Results.CVSubsegmentName;
            obj.CVSubsegmentNumber = parser.Results.CVSubsegmentNumber;
            obj.evaluationMetric = parser.Results.evaluationMetric;
            obj.intermediateResultFile = parser.Results.intermediateResultFile;
            obj.classificationPipeline = classificationPipeline;
            
            obj.parameterSearchFinished = false;
            
            
            
            obj.classNames = trainingPartition.classes.labels;
            
            obj.trainLabels = [];
            obj.testLabels = [];
            
            obj.prediction = [];
            obj.score = [];

            
            if ~iscell(parameters)
                parameters = {parameters};
            end
            
            obj.kFoldMatrix = MSKFoldMatrix(...
                trainingPartition.cvEntities.data);
            
            
            nParameters = length(parameters);
            nCVEntities = size(obj.kFoldMatrix,2);
            
            classNames = trainingPartition.classes.labels;

            % Create label vectors from MSLabelData.
            obj.trainLabels = cell(nCVEntities,1);
            obj.testLabels = cell(nCVEntities,1);
            for i=1:nCVEntities
                obj.trainLabels{i} = ...
                    trainingPartition.classes.data .* (obj.kFoldMatrix(:,i)==1);
                obj.testLabels{i} = ...
                    trainingPartition.classes.data .* (obj.kFoldMatrix(:,i)==2);
            end
            %% Create the different inner CV partitions.
            obj.CV = cell(nCVEntities,1);
            for i=1:nCVEntities
                
                obj.CV{i} = struct();
                obj.CV{i}.score = 0;
                obj.CV{i}.classNames = {};
                for j=1:length(classNames)
                    obj.CV{i}.classNames{j} = ...
                        ['CV_',num2str(i),'_T_',classNames{j}];
                end
                
                obj.CV{i}.classes = MSLabelData(...
                    obj.CV{i}.classNames,...
                    (obj.kFoldMatrix(:,i)==1).*trainingPartition.classes.data,...
                    obj.trainingData);
                
                obj.CV{i}.cvEntities = obj.CV{i}.classes.kFold(...
                    obj.CVSubsegmentNumber,...
                    trainingPartition.(obj.CVSubsegmentName));
                
                obj.CV{i}.kFoldMatrix = MSKFoldMatrix(...
                    obj.CV{i}.cvEntities.data);
                
                obj.CV{i}.trainLabels = cell(1,obj.CVSubsegmentNumber);
                obj.CV{i}.testLabels = cell(1,obj.CVSubsegmentNumber);
                for k=1:obj.CVSubsegmentNumber
                    obj.CV{i}.trainLabels{k} = ...
                        obj.CV{i}.classes.data .* (obj.CV{i}.kFoldMatrix(:,k)==1);
                    obj.CV{i}.testLabels{k} = ...
                        obj.CV{i}.classes.data .* (obj.CV{i}.kFoldMatrix(:,k)==2);
                end
                
                obj.CV{i}.P = cell(nParameters,1);
                
                obj.CV{i}.score = 0;
                obj.CV{i}.bestParameterIndex = 0;
                for j=1:nParameters
                    obj.CV{i}.P{j} = struct();
                    obj.CV{i}.P{j}.parameter = parameters{j};
                    for k=1:obj.CVSubsegmentNumber
                        obj.CV{i}.P{j}.S{k} = struct();
                    end
                end
                
            end
            
            
            
        end
        
        %% Aggregate the labels of a partitioned prediction
        function output = aggregateLabels(obj,results)
            aggregate = zeros(size(results{1}.prediction.data));
            for i=1:length(results)
                aggregate = ...
                    aggregate + results{i}.prediction.data;
            end
            output = MSLabelData(...
                obj.trainingPartition.classes.labels,...
                aggregate,...
                obj.trainingData);
        end
        
        %% Calculate the score based on evaluationMetric.
        function score = evaluate(obj,results)
            [~,confusionMatrix] = results.validate(...
                obj.trainingPartition.classes.data .* (results.data~=0));
            confusionMeasure = MSConfusionMeasure(confusionMatrix);
            score = confusionMeasure.(obj.evaluationMetric);
        end
        
        
        %% Conduct the whole study.
        function conductStudy(obj)
            if ~obj.parameterSearchFinished
                obj.parameterOptimization()
            end
            nCVEntities = length(obj.trainingPartition.cvEntities.labels);
            for i=1:nCVEntities
                obj.CV{i}.prediction = ...
                    obj.classificationPipeline(...
                    obj.trainingData,...
                    obj.trainLabels{i},...
                    obj.testLabels{i},...
                    obj.CV{i}.bestParameter,...
                    i,0);
                obj.CV{i}.predictionLabels = ...
                    obj.CV{i}.prediction.data;
                
            end
            obj.prediction = obj.aggregateLabels(obj.CV);
            obj.score = obj.evaluate(obj.prediction);
            disp(obj.score)
        end
        %% Parameter optimization.
        function parameterOptimization(obj)
            nCVEntities = length(obj.trainingPartition.cvEntities.labels);
            nParameters = length(obj.parameters);
            for i=1:nCVEntities
                for j=1:nParameters
                    for k=1:obj.CVSubsegmentNumber
                        obj.CV{i}.P{j}.S{k}.prediction = ...
                            obj.classificationPipeline(...
                            obj.trainingData,...
                            obj.CV{i}.trainLabels{k},...
                            obj.CV{i}.testLabels{k},...
                            obj.parameters{j},...
                            i,j);
                        obj.CV{i}.P{j}.S{k}.predictionLabels = ...
                            obj.CV{i}.P{j}.S{k}.prediction.data;
                    end
                    obj.CV{i}.P{j}.prediction = ...
                        obj.aggregateLabels(obj.CV{i}.P{j}.S);
                    
                    obj.CV{i}.P{j}.score = obj.evaluate(...
                        obj.CV{i}.P{j}.prediction);
                    disp(obj.CV{i}.P{j}.score);
                    
                    if obj.CV{i}.P{j}.score > obj.CV{i}.score
                        
                        
                        
                        obj.CV{i}.score = obj.CV{i}.P{j}.score;
                        
                        modifier = 0;
                        if (j+obj.parameterChoiceModifier < nParameters)&&...
                           (j+obj.parameterChoiceModifier > 0)    
                            modifier = obj.parameterChoiceModifier;
                        end
                        obj.CV{i}.bestParameterIndex = j+modifier;
                        obj.CV{i}.bestParameter = obj.parameters{j+modifier};
                    end
                end
            end
            obj.parameterSearchFinished = true;
        end
    end
    
end
