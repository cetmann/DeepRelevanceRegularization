classdef MSClassificationModel < matlab.mixin.Copyable
  % Abstract base class for classification model classes
  % 
  % This class is an abstract base class for classification models that map
  % data items (typically feature vectors) to label vectors. A derived
  % classification model has to implement the classify_impl() method that
  % computes the label vector for a set of data items. Moreover, a derived 
  % classification model shall define a unique creator string that will be 
  % stored in output label data as an indication of how the respective
  % labels were generated.
  %
  % Properties:
  %   creator:    Classification model creator
  %   classes:    Names of classes for which model was created
  %   numClasses: Number of classes
  %   modelInfo:  Struct for arbitrary classification model information
  %
  % Methods:
  %   MSClassificationModel: Constructor
  %   setClasses: Set class names
  %   classify: Apply classification model to input data
  %   classify_impl (abstract): Implementation of classification model
  %
  % MSClassificationModel uses the handle semantic, i.e. when assigning an
  % object of this class to a variable, only a reference to the original
  % object is copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = immutable)
    creator % Classification model creator

  end
  properties (SetAccess = private)
      classes % Names of classes for which model was created
  end

  properties (Dependent)
    numClasses % Number of classes
  end
  
  properties
    modelInfo = struct; % Struct for arbitrary classification model information
  end

  methods
    function obj = MSClassificationModel (creator)
      % Constructor
      % obj = MSClassificationModel(creator)
      %   creator: Optional model creator string
      
      % Check input arguments
      if nargin < 1
        creator = '';
      end
      if ~ischar(creator)
        error('creator must be either empty or a string')
      end
      obj.creator = creator;
    end
    
    function setClasses (obj, classes)
      % Set class names
      % obj.setClasses(classes): Set names of classes for which model was
      %   trained. classes must be a cell string vector.
      if ~(iscellstr(classes) && ~isempty(classes) && isvector(classes))
        error('classes must be a cell string vector')
      end
      obj.classes = classes;
    end
    function value=get.numClasses(obj)
        value=numel(obj.classes);
    end
  end
      
  methods
    function [outLabels, outScores] = classify (obj, msData, itemMask)
      % Apply classification model to input data
      % outLabels = obj.classify (msData, itemMask): 
      %   msData: Input data items to classify (typically feature data)
      %   itemMask: Optional logical or numerical vector or MSLabelData 
      %             object specifying a subset of data items as input data. 
      %             If specified, length must match msData.numItems.
      %   outLabels: MSLabelData object representing the predicted labels
      %              obtained from the classification model
      
      % Check input arguments
      narginchk(2,3)
      if ~isa(msData, 'MSData')
        error('Data argument must be an MSData object');
      end
      if nargin < 3
        % No annotations specified, select all data items
        itemMask = true(msData.numItems,1);
      elseif isa(itemMask, 'MSLabelData')
        if ~(itemMask.dataLength == 1 && itemMask.numItems == msData.numItems)
          error(['itemMask must have a single column, length must match msData.numItems ', ...
                 '(= %d)'], msData.numItems);
        else
          itemMask = itemMask.data > 0;
        end
      else
        if isnumeric(itemMask)
          itemMask = logical(itemMask);
        end
        if ~(islogical(itemMask) && isvector(itemMask) && ...
             length(itemMask) == msData.numItems)
          error(['itemMask must be a logical vector of length numItems ', ...
                 '(= %d)'], msData.numItems);
        end        
      end
      msData.assert;
      
      msData = MSClassifier.extractAdditionalData( msData );
      % Call classification model implementation
      L = zeros(msData.numItems,1);
      if nargout==1
          L(itemMask) = obj.classify_impl(msData, itemMask);
      else          
          [L(itemMask), scores]=obj.classify_impl(msData, itemMask);
          %if the classification methods returns scores
          if ~isempty(scores)
              S=zeros(msData.numItems,obj.numClasses,'like',msData.data);
              S(itemMask,:)=scores;
              outScores=MSScoreData(obj.classes,S,msData);
          else
              outScores=[];
          end    
      end
      % Create output labels
      outLabels = MSLabelData(obj.classes, L, msData);      
    end
  end
  
  methods (Abstract, Access = protected)
    [prediction, scores] = classify_impl (obj, msData, itemMask)
    % Implementation of classification model
    % prediction = obj.classify_impl(msData, itemMask): Compute predicted
    %   labels for input data msData (MSData object). itemMask specifies
    %   a logical vector of length msData.numItems specifiying the data
    %   items for which a prediction shall be computed. The output argument
    %   prediction shall contain predicted labels for the items selected by
    %   itemMask, its length shall equal the number of true elements of
    %   itemMask.
  end
end

