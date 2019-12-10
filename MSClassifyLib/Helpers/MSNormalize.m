classdef MSNormalize < MSFeatureMap
%     Normalizes MSData so that the mean in each dimension is 0 and the
%     standard deviation is 1.
%     Properties:
%     meanVector: the array of means
%     stdVector:  the array of standard deviations
%     
%     Methods:
%     MSNormalize: Constructor
%     map: normalizes a sample according to the learned means and stds
    
  properties
    meanVector % array of means
    stdVector %array of standard deviations
  end
  
  methods
    function obj = MSNormalize (msData,itemMask,creator)
      narginchk(1,2)
      if ~isa(msData, 'MSData')
        error('msData must be an MSData object')
      end
      if nargin < 2
        % No item mask given, select all items
        itemMask = true(msData.numItems,1);
      elseif ~(islogical(itemMask) && isvector(itemMask) && ...
          length(itemMask) == msData.numItems)
        error(['itemMask must be a logical vector of length numItems' ...
          '(= %d)'], msData.numItems)
      else
        itemMask = itemMask(:);
      end
      creator = 'MSNormalize';
      obj@MSFeatureMap(creator);
      obj.meanVector = mean(msData.data(itemMask,:),1);
      obj.stdVector = std(msData.data(itemMask,:),0,1);
    end
    function featureData = map (obj, msData,itemMask)
      narginchk(2,3)
      if ~isa(msData, 'MSData')
        error('msData must be an MSData object')
      end
      if nargin < 3
        % No item mask given, select all items
        itemMask = true(msData.numItems,1);
      elseif ~(islogical(itemMask) && isvector(itemMask) && ...
          length(itemMask) == msData.numItems)
        error(['itemMask must be a logical vector of length numItems' ...
          '(= %d)'], msData.numItems)
      else
        itemMask = itemMask(:);
      end
      data = bsxfun(@minus,msData.data(itemMask,:),obj.meanVector);
      data = bsxfun(@rdivide,data,obj.stdVector);
      data(find(isnan(data))) = 0;
      featureData = MSFeatureData(data, obj.creator, msData);   
    end
  end
end