classdef MSLDAClassifier < MSClassifier
  % Linear discrimant classifier (LDA)
  %
  % This class implements a linear discriminant classifier (LDA) for
  % multiple training classes. The generated classification model is an
  % MSLDAModel object.
  %
  % Methods:
  %   MSLDAClassifier: Constructor
  %   trainModel_impl: Implementation of LDA training algorithm
  %
  % MSLDAClassifier uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties
    ldaParams;
  end
    
  methods
    function obj = MSLDAClassifier (varargin)
      % Constructor
      % obj = MSLDAClassifier: Create new MSLDAClassifier object
      % obj = MSLDAClassifier(Name, Value, ...): Specify optional
      %   name-value-pair arguments that are passed to fitcdiscr.
      % Possible parameters are:
      % 'Prior': -string 'empirical' (according to the frequency of classes)
      %           or 'uniform'
      
       obj@MSClassifier;
       params=inputParser;

       params.addParameter('Prior', 'uniform', @checkPriorParam);
       listDT = {'linear','diaglinear','pseudolinear','quadratic','diagquadratic','pseudoquadratic'};
       params.addParameter('DiscrimType','linear', @(x)any(strcmp(x,listDT)));

       params.parse(varargin{:});

       obj.ldaParams = params.Results;
    end
  end
  
  methods (Access = protected)
    function model = trainModel_impl (obj, msData, L, numFeatures)
      % model = obj.trainModel_impl(msData, L, numFeatures): Create 
      %   classification model for input data msData (MSData) and classes
      %   defined by the labels argument (numeric vector). Only the first
      %   numFeatures components (columns) of the input data are used.

      % Train LDA
      nonZeroLabels = (L > 0);
      ldaParamsNVPair = struct2NameValuePair(obj.ldaParams);
      lda = fitcdiscr( msData.data(nonZeroLabels, 1:numFeatures), ...
              L(nonZeroLabels), ldaParamsNVPair{:} );
      % Remove training data from LDA model
      lda = lda.compact();
      
      model = MSLDAModel(lda, numFeatures, class(obj));
    end
  end
  
end
function bool = checkPriorParam(param)
bool = ischar(param)||(isnumeric(param) && isvector(param))||...
       (isstruct(param) && isfield(param, 'ClassProbs') && isfield(param, 'ClassNames') );
end
