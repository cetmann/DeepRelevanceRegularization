classdef MSLDAModel < MSClassificationModel
  % Linear discriminant analysis (LDA) classification model
  %
  % This class implements the linear discriminant analysis (LDA)
  % classification model..
  %
  % Properties:
  %   lda: LDA model as created by fitcdiscr()
  %   numFeatures: Number of features used for training
  %
  % Methods:
  %   MSLDAModel: Constructor
  %   classify_impl: Classification model implementation
  %
  % MSLDAModel uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = immutable)
    lda;         % LDA model as created by fitcdiscr()
    numFeatures; % Number of features used for training
  end
  
  methods
    function obj = MSLDAModel (lda, numFeatures, creator)
      % Constructor
      obj@MSClassificationModel(creator);
      obj.lda = lda;
      obj.numFeatures = numFeatures;
    end
  end
  
  methods (Access = protected)
    function [prediction, scores] = classify_impl (obj, msData, itemMask)
      % Classification model implementation
      
      % Compute prediction for items in itemMask
      if nargout>1
          [prediction, scores] = obj.lda.predict(msData.data(itemMask,1:obj.numFeatures));
      else
          prediction = obj.lda.predict(msData.data(itemMask,1:obj.numFeatures));
      end
    end    
  end
  
end

