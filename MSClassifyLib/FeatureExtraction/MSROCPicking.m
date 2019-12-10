classdef MSROCPicking < MSFeatureExtraction
  % Feature extraction based on supervised mz-values selection by ROC criterion
  %
  % This class implements a supervised feature extraction method based on
  % computing ROC values (area under receiver operating characteristics
  % curve) for each mz-value. For two classes, the ROC value represents the
  % probability that values in one class are larger than values in the
  % other class. For more than two classes, the maximum over multiple 
  % one-versus-rest ROC values is computed.
  %
  % Properties:
  %    rankCrit: Criteria used to select the most relevant peaks
  %
  % Methods:
  %   MSROCPicking: Constructor
  %   createMap_impl: Create mz-values map from ROC values
  %
  % MSROCPicking uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = protected)
    rankCrit  % 
  end
  
  methods
    function obj = MSROCPicking (maxNumFeatures, rankCrit)
      % Constructor
      % obj = MSROCPicking(maxNumFeatures):      
      %   maxNumFeatures:  Maximum overall number of features used in map
      
      % Check input arguments
      narginchk(1,2);
      if nargin<2
        rankCrit = 'rocMax';
      elseif ~any(strcmp(rankCrit,{'rocMax','rocMaxLabelBalance'}))
        error('Unknown rank criterion');
      end
      if nargin < 1
        maxNumFeatures = inf;
      end
      % Initialize properties
      obj@MSFeatureExtraction(maxNumFeatures);
      obj.rankCrit=rankCrit;
    end
  end
  methods (Access=protected)
    function map = createMap_impl (obj, msData, labels)
      % Create mz-values map from ROC values with respect to class labels
      % map = obj.createMap (msData, labels): Create mz-values map for
      %   spectral dataset msData (MSData) and classes defined by the
      %   labels argument (numeric labels vector)
      %
      %   The created map object is an MSMzValuesMap, the featureMapInfo
      %   property contains a complete list of the ranked ROC values.

      
      % get ranked list and rocValues matrix (to be more exact structs 
      % with .rocMax' and .rocMaxLabelBalance fields) 
      if length(unique(labels(labels~=0)))< 2
          error('Not enough labels specified to create the feature map')
      end
      [rankedList,rocMaxVals,rocValsTotal,pValsMax, pValsTotal, discrimClass] = ...
                                        MSROCRanking(msData,labels,'true');
      
      % Create mz-values feature map
      indexList = rankedList.(obj.rankCrit);
      map = MSMzValuesMap(msData.mzVector, indexList(1:obj.maxNumFeatures), ...
                          class(obj));                        
      % Store full list of ROC values as additional feature map information
      map.featureMapInfo.rankCrit = obj.rankCrit;
      map.featureMapInfo.rocMaxVals = rocMaxVals;
      map.featureMapInfo.rocValsTotal = rocValsTotal;
      map.featureMapInfo.rankedList = rankedList;
      map.featureMapInfo.pValsMax = pValsMax;
      map.featureMapInfo.pValsTotal = pValsTotal;
      map.featureMapInfo.discrimClass = discrimClass;
    end
  end
end





