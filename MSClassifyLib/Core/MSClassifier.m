classdef MSClassifier < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
  % Abstract base class for classifier classes
  % 
  % This class is an abstract base class for classifier classes that create
  % a classification model trained for a given set of data and training 
  % classes. A derived classifier has to define the trainModel_impl() 
  % method that creates a classification model object.
  %
  % Properties:
  %   supportScores: bool specifying whether the classifier supports score
  %     classification
  %
  % Methods:
  %   trainModel: Perform classifier training on input data
  %   trainModel_impl (abstract): Implementation of training algorithm
  %
  % MSClassifier uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = protected)
      supportScores = true;
  end
  
  methods
    function model = trainModel (obj, msData, labels, numFeatures)
      % Perform classifier training on input data
      % model = obj.trainModel(msData, labels): Create classification model
      %   for input data msData (MSData) and classes defined by the labels 
      %   argument (MSLabelData or numeric labels vector). The optional 
      %   argument numFeatures may be specified to limit the number of data
      %   components (columns) used for training the classification model.
      
      % Check input data argument
      if ~isa(msData, 'MSData')
        error('Input data argument must be an MSData object');
      end
      msData.assert;
      if nargin < 4
        numFeatures = [];
      end
      if ~isempty(numFeatures) && ~isscalar(numFeatures)
        error('Argument numFeatures must be a positive integer scalar');
      end
      if isempty(numFeatures)
        % By default use all available features
        numFeatures = msData.dataLength;
      else
        % Number of features limited by available features
        numFeatures = min(floor(numFeatures), msData.dataLength);
      end
      
      % Check labels argument and create label vector
      L = MSLabelData.getLabelVector(labels, msData);
      % Call classifier implementation
      model = obj.trainModel_impl(obj.extractAdditionalData(msData), L, numFeatures);
      % Set class names
      if isa(labels, 'MSLabelData')
        model.setClasses(labels.labels);
      else
        maxLabel = floor(max(labels));
        if maxLabel >= 1
          classes = cellstr(num2str((1:maxLabel)','%-g'));
        else
          classes = {};
        end
        model.setClasses(classes);
      end
    end
  end  
  
  methods (Abstract, Access = protected)
    model = trainModel_impl (obj, msData, L, numFeatures)
    % Implementation of training algorithm
    % model = obj.trainModel_impl(msData, L, numFeatures): Create 
    %   classification model for input data msData (MSData) and classes
    %   defined by the labels argument (numeric vector). Only the first
    %   numFeatures components (columns) of the input data are used.
  end
  
  methods (Static)
      function newMsData = extractAdditionalData( msData )
         nAddD = length( msData.additionalDataVersions );
         if nAddD > 0
             newMsData = msData.copy;
             for i = 1:nAddD
                newMsData.data = [ newMsData.data, msData.additionalDataVersions{i}.data ]; 
             end
         else
             newMsData = msData;
         end
      end
  end
  
end

