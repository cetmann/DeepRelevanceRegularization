classdef MSMzValuesMap < MSFeatureMap
  % Feature map based on picking specific data values
  % 
  % This class is a simple feature map picking input data values at
  % specific index positions and copying these to the output feature
  % vectors. An MSmzValuesMap is parameterized by index vector and the
  % expected data length of an input data object. The mz-values expected at
  % the index positions may also be specified. An error is generated in
  % case the input data properties do not match.
  %
  % Properties:
  %   index:       Vector of indices used for picking output values
  %
  % Methods:
  %   MSMzValuesMap: Constructor
  %   map_impl: Compute output feature data from input data
  %
  % MSmzValuesMap uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = protected)
    index      % Vector of indices used for picking output values
  end
  
  methods
    function obj = MSMzValuesMap (mzVector, index, creator)
      % Constructor
      % obj = MSMzValuesMap(mzVector, index, creator)
      %   mzVector: m/z vector of expected input data items
      %   index: Vector of indices used for picking output values
      %   creator: Feature map creator string (optional, defaults to MSmzValuesMap)
      
      narginchk(2,3);
      if nargin < 3
        % If creator string is not specified, use this class name
        creator = 'MSMzValuesMap';
      end
      
      % Check input arguments
      if ~(~isempty(index) && isnumeric(index) && isvector(index))
          error('index must be an integer vector');
      end
      % Call super class constructor where mzVector validity is checked
      obj@MSFeatureMap(mzVector, creator);
      % Check whether specified index values are unique and in proper range
      u = unique(index);
      if ~(u(1) >= 1 && u(end) <= length(obj.mzVector) && length(u) == length(index))
        error('index must contain unique index values in [1..%d]',length(mzVector));
      end         
      obj.index = floor(index);
    end 
    
    function [auc, pvalue, discrimClass] = maxAUC(obj)
        info = obj.featureMapInfo;
        auc = info.rocMaxVals(obj.index);
        pvalue = info.pValsMax(obj.index);
        discrimClass = info.discrimClass(obj.index);
    end
    
    function show(obj, type, varargin)
        narginchk(1,inf)
        typeList = {'auc','pvalue'};
        if nargin < 2
            type = 'auc';
        elseif ~any(strcmp(type,typeList))
            error('Unknown type of graphic')
        end
        [aucValues, pValues, discrimClass] = obj.maxAUC;
        classes = unique(discrimClass);
        if strcmp(type, 'auc')
            bgValues = obj.featureMapInfo.rocMaxVals;
            fgValues = aucValues;
        else
            bgValues = obj.featureMapInfo.pValsMax;
            fgValues = pValues;
        end
        % plot in gray all roc values
        gray = [211/255,211/255,211/255];
        bar(obj.mzVector,bgValues,'EdgeColor',gray,...
            'FaceColor', gray), hold on, axis tight
        nclasses = length(classes);
        v = zeros(1,nclasses);
        % plot best roc values (highest absolute value)
        for i=1:nclasses
            mask = discrimClass==classes(i);
            v(i) = stem(obj.mzVector(obj.index(mask)),...
                fgValues(mask),'filled','linewidth', 1); hold on
        end
        hold off
        if nclasses > 1|| any(classes > 1)
            legend(v, num2str(classes));
        end
    end
  end
    
  methods (Access=protected)
    function featureData = map_impl (obj, msData, itemMask, numFeatures, ~)
      % Compute output feature data from input data
      % featureData = obj.map(msData, itemMask, numFeatures, printFlag)
      %   msData: MSData object with input data (spectra, feature vectors)
      %   featureData: Resulting MSFeatureData object with feature vectors
      %                consisting of input data values at the index
      %                positions defined for the feature map.
      %   itemMask: logical mask specifying the spectra to which the mapping
      %             is applied
      %   numFeatures: Restrict the output number of features to the
      %                specified quantity
      %   printFlag: ignored
      
      % Create empty output feature data
      featureData = MSFeatureData(zeros(msData.numItems, numFeatures, 'like', msData.data), ...
                                  class(obj), msData, itemMask);
      % Compute feature data for items selected by itemMask
      featureData.data(itemMask,:) = msData.data(itemMask,obj.index(1:numFeatures));
      % Store mz-values from input object into featureInfo of output data
      featureData.featureInfo.mzValues = msData.mzVector(obj.index(1:numFeatures));
    end
    
    function numFeatures = getNumFeatures(obj)
      % Get number of output features  
      numFeatures = length(obj.index);
    end    
  end
end

