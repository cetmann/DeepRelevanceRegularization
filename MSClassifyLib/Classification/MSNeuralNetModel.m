classdef MSNeuralNetModel < MSClassificationModel
  % Neural network classification model
  properties (SetAccess = immutable)
    
    numFeatures 
    % Number of used features.
    
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
    
    resultFile
    % File used for data exchange between Python and
    % Matlab. Contains predicted class labels and
    % predicted class membership probabilites.
    
    pyVersion
    % Version of your python installation
    % You can choose from 'python', 'python2' or 'python3'
    
    backend
    % The backend that is used for the experiment
    % One can choose between the classic 'theano' and the
    % more modern 'keras' backend.
    % Note that 'theano' might not be supported in the future
    % due to end-of-life of most packages


  end
  
  properties
    inputDimScaling
    % scling of input dimension via standard deviation of each input
    % dimension of the used training data. Scaling parameter is used to
    % scale sensitivity output   
    
    keepTemporaryFiles
    % true if temporary files should be kept, else false
  end
  
  methods
    function obj = MSNeuralNetModel (numFeatures,...
        inputDimScaling,...
        neuralNetFile,neuralNetArchitecture,...
        dataFile,resultFile,... 
        pyVersion, backend,...
        keepTemporaryFiles, creator)
      % Constructor
      obj@MSClassificationModel(creator);
      obj.numFeatures = numFeatures;
      obj.inputDimScaling = inputDimScaling;
      obj.neuralNetFile = neuralNetFile;
      obj.neuralNetArchitecture = neuralNetArchitecture;
      obj.dataFile = dataFile;
      obj.resultFile = resultFile;
      obj.pyVersion = pyVersion;
      obj.backend = backend;
      obj.keepTemporaryFiles = keepTemporaryFiles;
    end
    
  end
  
  methods (Access = protected)
    function [prediction, scores] = classify_impl (obj, msData, itemMask)
      % prepare data file which is load by python scripts
      prepareDataFile(obj, msData, itemMask, false)

      % The following lines call a python program which performs the
      % network's predicting.
      %pythonPrefix = '~/NeuralNets/env/bin/';
      pythonPrefix = {''};
      pyCommand = strcat(pythonPrefix,...
        obj.pyVersion,...
        {' -B '},...
        which('NNInterface.py'),...
        {' -neuralNetFile '}, ...
        obj.neuralNetFile,...
        {' -dataFile '}, ...
        obj.dataFile,...
        {' -resultFile '}, ...
        obj.resultFile,...
        {' -backend '}, ...
        obj.backend, ...
        {' -mode '},{'predict'});
      
      system(pyCommand{1});
      
      % Compute prediction for items in itemMask
      
      if nargout > 1
        scores = h5read(obj.resultFile,'/scores');
      end
      prediction = h5read(obj.resultFile,'/prediction');
      

      %[~,prediction]=max(prediction_temp.probabilities,[],2);
      if ~obj.keepTemporaryFiles
        delete(obj.dataFile,...
          obj.resultFile)
      end
    end
    
    function [data] = prepareDataFile(obj, msData, itemMask, transferLabels)
      if isempty(obj.numFeatures)
        N = msData.dataLength;
      else
        N = obj.numFeatures;
      end
      
      % get mask if itemMask was provided as label data
      if isa(itemMask, 'MSLabelData')
          itemMask = itemMask.data;
      end
      
      % reduce data to mask
      data = msData.data(itemMask>0,1:N);
      data(isnan(msData.data)) = 0;
      if transferLabels
        % add labels if used for NLL sensitivty analysis
          classes = itemMask(itemMask > 0, 1);
          save(obj.dataFile, 'data', 'classes','-v7.3');
          clearvars data classes
      else 
          save(obj.dataFile, 'data','-v7.3');
          clearvars data 
      end

    end
    function [] = prepareDataFileReconBasis(obj, msData, itemMask, basisRecon)
      if isempty(obj.numFeatures)
        N = msData.dataLength;
      else
        N = obj.numFeatures;
      end
      
      % get mask if itemMask was provided as label data
      if isa(itemMask, 'MSLabelData')
          itemMask = itemMask.data;
      end
      
      % reduce data to mask
      data = msData.data(itemMask>0,1:N);
      data(isnan(msData.data)) = 0;
      disp(size(data))
      save(obj.dataFile, 'data', 'basisRecon', '-v7.3');
      clearvars data basis

    end
    
  end
  
  methods (Access = public)
    function [sensitivity] = sensitivityAnalysis(obj, msdata, itemMask, type)
        
%             parser = inputParser;
% 
%             @isstr does not work for input parsers.
%             checkString = @(x) validateattributes(x,{'char'},{'nonempty'});
%             addParameter(parser,'trainingDataStd','data_std',checkString);
%             
%             parse(parser,varargin{:});
%             
%             obj.hyperparameters = struct();
%             obj.hyperparameters.trainingDataStd = ...
%                 parser.Results.trainingDataStd;
      
      % select type of sensitivity and prepare data file
      if strcmp(type, 'NLL') % sensitivity of nonnegative log-likelihood
          prepareDataFile(obj, msdata, itemMask, true)
      else if strcmp(type, 'perClass') || strcmp(type, 'perClassLogit')% sensitivity per class
              prepareDataFile(obj, msdata, itemMask, false)
          else % sensitivity of binary crossentropy  
              type = 'binary';
              prepareDataFile(obj, msdata, itemMask, false)
          end
      end
      
      % The following lines call a python program which performs the
      % network's sensitivity analysis.
      pythonPrefix = '';
      pyCommand = strcat(pythonPrefix,...
        obj.pyVersion,...
        {' -B '},...
        which('NNInterface.py'),...
        {' -neuralNetFile '}, ...
        obj.neuralNetFile,...
        {' -dataFile '}, ...
        obj.dataFile,...
        {' -resultFile '}, ...
        obj.resultFile,...
        {' -backend '}, ...
        obj.backend, ...
        {' -mode '}, {'sensitivity'}, ...
        {' -sensType '}, {type}); % Compute sensitivity for items in itemMask
      
      
      % call python code from console
      system(pyCommand{1});
      
      % Load result file with sensitivity computations
      sensitivity = squeeze(h5read(obj.resultFile,'/sensitivity'));
      
      % rearange dimension
      if strcmp(type,'perClass') || strcmp(type, 'perClassLogit')
        sensitivity = permute(sensitivity, [3,2,1]);
      end
      
      % scale sensitivity by property inputDimScaling
      if length(size(sensitivity)) > 2 % if sensitivity per class
          disp("!!!!!")
          for i = 1: size(sensitivity, 1)
              sensitivity(i, :, :) = bsxfun(@times,squeeze(sensitivity(i, :, :)), obj.inputDimScaling);
          end
      else
        sensitivity = bsxfun(@times,sensitivity, obj.inputDimScaling);
      end
      
      % delete temporary result file if desired
      if ~obj.keepTemporaryFiles
        delete(obj.dataFile,...
          obj.resultFile)
      end  
    end
    
    function [reconInput] = reconstructInput(obj, msdata, itemMask, varargin)
        % parse varargin
        parser = inputParser;
            
        isPositiveNumeric = @(x) isnumeric(x) && x>0;
        isNonnegativeNumeric = @(x) isnumeric(x) && x>=0;
        
        addParameter(parser,'weighting',1,@isnumeric);
        addParameter(parser,'layer',1,isNonnegativeNumeric);
        addParameter(parser,'stepSize',0.1,isPositiveNumeric);
        addParameter(parser,'numIterations',10,isPositiveNumeric);
        addParameter(parser,'predCoeff',0.1,isNonnegativeNumeric);
        addParameter(parser,'regCoeffL1',0.1,isNonnegativeNumeric);
        addParameter(parser,'regCoeffL2',0.1,isNonnegativeNumeric);            
            
        parse(parser,varargin{:});
        % create object with hyperparameter
        hyperparameters = struct();
        hyperparameters.weighting = ...
            parser.Results.weighting;
        hyperparameters.reconLayer = ...
            int32(ceil(parser.Results.layer));
        hyperparameters.reconStepSize = ...
            parser.Results.stepSize;
        hyperparameters.reconNumIterations = ...
            int32(ceil(parser.Results.numIterations));
        hyperparameters.predCoeff = ...
            parser.Results.predCoeff;
        hyperparameters.reconRegCoeffL1 = ...
            parser.Results.regCoeffL1;
        hyperparameters.reconRegCoeffL2 = ...
            parser.Results.regCoeffL2;          
            
      % create data file
      prepareDataFile(obj, msdata, itemMask, true)
      
      % The following lines call a python program which performs the
      % input reconstruction for given data.
      pythonPrefix = '';
      pyCommand = strcat(pythonPrefix,...
        obj.pyVersion,...
        {' -B '},...
        which('NNInterface.py'),...
        {' -neuralNetFile '}, ...
        obj.neuralNetFile,...
        {' -dataFile '}, ...
        obj.dataFile,...
        {' -resultFile '}, ...
        obj.resultFile,...
        {' -backend '}, ...
        obj.backend, ...
        {' -mode '}, {'reconstruction'});
    
      % Add hyperparameter to command line
      hyperparameterKeys = fieldnames(hyperparameters)';
            
      for key = hyperparameterKeys
          pyCommand = strcat(pyCommand, {[' -',char(key)]});
          pyCommand = strcat(pyCommand, ...
             {[' ',num2str(getfield(hyperparameters,char(key)))]});
      end
      
      % call python code from console
      system(pyCommand{1});
      
      % Load result file with sensitivity computations
      reconInput = squeeze(h5read(obj.resultFile,'/reconstruction'));
      
      % delete temporary result file if desired
      if ~obj.keepTemporaryFiles
        delete(obj.dataFile,...
          obj.resultFile)
      end  
    end
    
    function [reconInput] = reconstructInputinBasis(obj, msdata, itemMask, ...
                                                    basis, varargin)
        % parse varargin
        parser = inputParser;
            
        isPositiveNumeric = @(x) isnumeric(x) && x>0;
        isNonnegativeNumeric = @(x) isnumeric(x) && x>=0;
        
        addParameter(parser,'layer',1,isNonnegativeNumeric);
        addParameter(parser,'stepSize',0.1,isPositiveNumeric);
        addParameter(parser,'numIterations',10,isPositiveNumeric);
        addParameter(parser,'regCoeffL1',0.1,isNonnegativeNumeric);           
            
        parse(parser,varargin{:});
        % create object with hyperparameter
        hyperparameters = struct();
        hyperparameters.reconLayer = ...
            int32(ceil(parser.Results.layer));
        hyperparameters.reconStepSize = ...
            parser.Results.stepSize;
        hyperparameters.reconNumIterations = ...
            int32(ceil(parser.Results.numIterations));
        hyperparameters.reconRegCoeffL1 = ...
            parser.Results.regCoeffL1;
            
      % create data file
      prepareDataFileReconBasis(obj, msdata, itemMask, basis)
      
      % The following lines call a python program which performs the
      % input reconstruction for given data.
      pythonPrefix = '';
      pyCommand = strcat(pythonPrefix,...
        obj.pyVersion,...
        {' -B '},...
        which('NNInterface.py'),...
        {' -neuralNetFile '}, ...
        obj.neuralNetFile,...
        {' -dataFile '}, ...
        obj.dataFile,...
        {' -resultFile '}, ...
        obj.resultFile,...
        {' -backend '}, ...
        obj.backend, ...
        {' -mode '}, {'reconstructionBasis'});
    
      % Add hyperparameter to command line
      hyperparameterKeys = fieldnames(hyperparameters)';
            
      for key = hyperparameterKeys
          pyCommand = strcat(pyCommand, {[' -',char(key)]});
          pyCommand = strcat(pyCommand, ...
             {[' ',num2str(getfield(hyperparameters,char(key)))]});
      end
      
      % call python code from console
      system(pyCommand{1});
      
      % Load result file with sensitivity computations
      reconInput = squeeze(h5read(obj.resultFile,'/reconstructionBasis'));
      
      % delete temporary result file if desired
      if ~obj.keepTemporaryFiles
        delete(obj.dataFile,...
          obj.resultFile)
      end  
    end
  end
  
end