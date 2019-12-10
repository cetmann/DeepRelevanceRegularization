classdef MSNeuralNetClassifier < MSClassifier
    properties
        %numFeatures
        % Number of features to be used for training
        
        pyVersion
        % Version of your python installation
        % You can choose from 'python', 'python2' or 'python3'
        % The standard value is 'python3'
        
        neuralNetFile
        % Serialized Neural Net. Also contains class
        % dictionaries that can translate between the class
        % labels used by Python (0,..,N-1) and the class
        % labels transmitted by MSClassifyLib (e.g. 1,3 and 4).
        
        neuralNetArchitecture
        % An external network architecture can be supplied.
        % If neuralNetArchitecture is 'auto', the standard architecture
        % defined in NNInterface.py is used.
        
        dataFile
        % File used for data exchange between Python and
        % Matlab. Contains spectra (and class labels during
        % training).
        % Default: 'neuralNetTemp'
        
        resultFile
        % File used for data exchange between Python and
        % Matlab. Contains predicted class labels and
        % predicted class membership probabilites.
        % Default: 'neuralNetTemp'
        
        backend
        % The backend that is used for the experiment
        % One can choose between the classic 'theano' and the
        % more modern 'keras' backend.
        % Note that 'theano' might not be supported in the future
        % due to end-of-life of most packages
        % Default: 'theano'
        
        continueTraining
        % false if the neuralNetFile should be created anew
        % every time the neural network is trained,
        % otherwise true.
        % Default: false
        
        keepTemporaryFiles
        % true if temporary files should be kept, else false
        
        useValidationSet
        % If true, 10% of the training data is used for validation.
        
        hyperparameters
        % Struct of hyperparameters.
        
        useTimeStamp
        % If 1, a timestamp is added to the neuralNetFile-name upon the
        % beginning of the training. If for example an 
        % MSClassificationStudy is performed with cross-validation, the
        % file name would be identical for all cross-validation folds, so
        % that each neuralNetFile would be overwritten.
        
        saveBestValidModel
        % If 1, the best model determined by the lowest valid loss during
        % training is saved and the final returned model is this best model
        %  default: 1
        
        sensScaling
        % scaling of input dimension via standard deviation of each input
        % dimension of the used training data. Scaling parameter is used to
        % scale sensitivity output
        
        useSensRegControl
        
        weightClasses
        % Whether to weight the classes according to their prevalence.
    end
    
    methods
        function obj = MSNeuralNetClassifier (varargin)
            
            obj@MSClassifier;
            parser = inputParser;
            
            
            % @isstr does not work for input parsers.
            checkString = @(x) validateattributes(x,{'char'},{'nonempty'});
            
            addParameter(parser,'pyVersion','python3',checkString)
            addParameter(parser,'neuralNetFile','neuralNetFile.nnet',checkString);
            addParameter(parser,'dataFile','dataFile.mat',checkString);
            addParameter(parser,'resultFile','resultFile.mat',checkString);
            addParameter(parser,'neuralNetArchitecture','auto',checkString);
            addParameter(parser,'backend','theano',checkString);
            addParameter(parser,'continueTraining',0,@isnumeric);
            addParameter(parser,'keepTemporaryFiles',0,@isnumeric);
            addParameter(parser,'useValidationSet',1,@isnumeric);
            addParameter(parser,'useTimeStamp',0,@isnumeric);
            addParameter(parser,'useSensRegControl',0,@isnumeric);
            addParameter(parser,'saveBestValidModel',0,@isnumeric);
            addParameter(parser,'sensScaling',0,@isnumeric);
            
            isPositiveNumeric = @(x) isnumeric(x) && x>0;
            isNonnegativeNumeric = @(x) isnumeric(x) && x>=0;
            isRatio = @(x) isnumeric(x) && x>=0 && x<=1;
            isZeroOrOne = @(x) x==0 || x==1;
            
            addParameter(parser,'weightClasses',1,isZeroOrOne);
            addParameter(parser,'verbose',2,isNonnegativeNumeric);
            addParameter(parser,'batchSize',128,isPositiveNumeric);
            addParameter(parser,'epochs',100,isPositiveNumeric);
            addParameter(parser,'lrSchedule',0,isNonnegativeNumeric);
            addParameter(parser,'sdWidth',0,isNonnegativeNumeric);
            addParameter(parser,'sdNumber',0,isNonnegativeNumeric);
            addParameter(parser,'learningRate',0.001,isPositiveNumeric);
            addParameter(parser,'adversarialNoise',0.,isNonnegativeNumeric);
            addParameter(parser,'l1',0.,isNonnegativeNumeric);
            addParameter(parser,'l2',0.,isNonnegativeNumeric);
            addParameter(parser,'logitSens',0.,isNonnegativeNumeric);
            addParameter(parser,'logitDiffSens',0.,isNonnegativeNumeric);
            addParameter(parser,'logitSqSens',0.,@isnumeric);
            addParameter(parser,'probSens',0.,isNonnegativeNumeric);
            addParameter(parser,'lossSens',0.,isNonnegativeNumeric);
            addParameter(parser,'disturbLabelRate',0,isRatio);
            addParameter(parser,'validationSetRatio',0.1,isRatio);
            addParameter(parser,'trainingDataStd','data_std',checkString);
            addParameter(parser,'lossFunc','categorical_crossentropy',checkString);
            
            parse(parser,varargin{:});
            
            obj.weightClasses = (parser.Results.weightClasses~=0);
            obj.pyVersion = parser.Results.pyVersion;
            obj.neuralNetFile = parser.Results.neuralNetFile;
            obj.dataFile = parser.Results.dataFile;
            obj.neuralNetArchitecture = parser.Results.neuralNetArchitecture;
            obj.backend = parser.Results.backend;
            obj.resultFile = parser.Results.resultFile;
            obj.continueTraining = (parser.Results.continueTraining~=0);
            obj.keepTemporaryFiles = (parser.Results.keepTemporaryFiles~=0);
            obj.useValidationSet = (parser.Results.useValidationSet~=0);
            obj.useTimeStamp = (parser.Results.useTimeStamp~=0);
            obj.useSensRegControl = (parser.Results.useSensRegControl~=0);
            obj.saveBestValidModel = (parser.Results.saveBestValidModel~=0);
            obj.sensScaling = parser.Results.sensScaling;
            
            obj.hyperparameters = struct();
            obj.hyperparameters.verbose = parser.Results.verbose;
            obj.hyperparameters.batchSize = ...
                int32(ceil(parser.Results.batchSize));
            obj.hyperparameters.epochs = ...
                int32(ceil(parser.Results.epochs));
            obj.hyperparameters.sdWidth = ...
                int32(ceil(parser.Results.sdWidth));
            obj.hyperparameters.sdNumber = ...
                int32(ceil(parser.Results.sdNumber));
            obj.hyperparameters.learningRate = ...
                parser.Results.learningRate;
            obj.hyperparameters.lrSchedule = ...
                parser.Results.lrSchedule;
            obj.hyperparameters.adversarialNoise = ...
                parser.Results.adversarialNoise;
            obj.hyperparameters.l1 = ...
                parser.Results.l1;
            obj.hyperparameters.l2 = ...
                parser.Results.l2;
            obj.hyperparameters.logitSens = ...
                parser.Results.logitSens;
            obj.hyperparameters.logitDiffSens = ...
                parser.Results.logitDiffSens;
            obj.hyperparameters.logitSqSens = ...
                parser.Results.logitSqSens;
            obj.hyperparameters.probSens = ...
                parser.Results.probSens;
            obj.hyperparameters.lossSens = ...
                parser.Results.lossSens;
            obj.hyperparameters.disturbLabelRate = ...
                parser.Results.disturbLabelRate;
            obj.hyperparameters.validationSetRatio = ...
                parser.Results.validationSetRatio;
            obj.hyperparameters.trainingDataStd = ...
                parser.Results.trainingDataStd;
            obj.hyperparameters.lossFunc = ...
                parser.Results.lossFunc;
        end
        
    end
    
    methods (Access = protected)
        function model = trainModel_impl (obj, msData, labels, numFeatures)
            % Number of features defaults to length of data items
            N = numFeatures;
            if isempty(N)
                N = msData.dataLength;
            end
            
            if N < msData.dataLength
                disp(['Used only the ',num2str(N),...
                    ' first features for training!'])
            end
            nonZeroLabels = (labels > 0);
            
            % If useTimeStamp is true, add a timestamp to the file name.
            if obj.useTimeStamp
               [pathstr,name,ext] = fileparts(obj.neuralNetFile);
               [name,~] = strsplit(name,'@'); 
               name = name{1};
               timestamp = ...
                   char(datetime('now','Format','yyyyMMdd_HHmmss'));
               obj.neuralNetFile = ...
                   fullfile(pathstr,[strcat(name,'@',timestamp) ext]);
            end
            
            % Use only labeled data for transfer
            data = msData.data(nonZeroLabels, 1:N);
            data(isnan(msData.data)) = 0;

            % Scale the saliency by a presupplied vector. If none is 
            % provided, use the data's standard deviation.
            if length(obj.sensScaling) ~= size(data,2)
                inputDimScaling = std(data,1);
                
                % If sensScaling is 1, produce a 1-vector of appropriate
                % size.
                if obj.sensScaling == 1
                    inputDimScaling = ones(size(inputDimScaling));
                    disp('Using no sensScaling.')
                % Only display this message is sparse saliency
                % regularization is actually applied.
                elseif obj.hyperparameters.logitSens | obj.hyperparameters.logitDiffSens | obj.hyperparameters.logitSqSens | obj.hyperparameters.probSens | obj.hyperparameters.lossSens
                   disp('No valid sensScaling provided. Using standard deviation.') 
                end
               
            end
            
            classes = labels(nonZeroLabels,1);
            
            save(obj.dataFile, 'data', 'classes','-v7.3');
            
            clearvars data classes;
            
            % Add custom arguments before the call of the python
            % environment
            pythonPrefix = {''};
            % -B prevents the creation of the bytecode-file (.pyc)
            pythonArguments = {' -B '};
            
            % The following lines call a python program which performs the
            % network's training.
            pyCommand = strcat(pythonPrefix,...
                obj.pyVersion,...
                pythonArguments,...
                which('NNInterface.py'),...
                {' -neuralNetFile '}, obj.neuralNetFile,...
                {' -dataFile '}, obj.dataFile,...
                {' -resultFile '}, obj.resultFile,...
                {' -backend '}, obj.backend,...
                {' -neuralNetArchitecture '}, obj.neuralNetArchitecture,...
                {' -mode '}, 'train');
            
            if ~obj.continueTraining
                pyCommand = strcat(pyCommand, {' -overwrite '});
            end
            
            if obj.useValidationSet
                pyCommand = strcat(pyCommand, {' -useValidationSet'});
            end
            
            if obj.saveBestValidModel
                pyCommand = strcat(pyCommand, {' -saveBestValidModel'});
            end
            
            if obj.useSensRegControl
                pyCommand = strcat(pyCommand, {' -useSensRegControl'});
            end
            if obj.weightClasses
                pyCommand = strcat(pyCommand, {' -weightClasses'});
            end
            
            hyperparameterKeys = fieldnames(obj.hyperparameters)';
            
            for key = hyperparameterKeys
                pyCommand = strcat(pyCommand, {[' -',char(key)]});
                pyCommand = strcat(pyCommand, ...
                    {[' ',num2str(getfield(obj.hyperparameters,char(key)))]});
            end
            
            system(pyCommand{1});
            
            model = MSNeuralNetModel(numFeatures,...
                inputDimScaling,...
                obj.neuralNetFile, obj.neuralNetArchitecture,...
                obj.dataFile, obj.resultFile,...
                obj.pyVersion, obj.backend,...
                obj.keepTemporaryFiles,class(obj));
            
            if ~obj.keepTemporaryFiles
                delete(obj.dataFile)
            end
        end
    end
end