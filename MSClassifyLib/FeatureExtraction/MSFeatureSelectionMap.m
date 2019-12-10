classdef MSFeatureSelectionMap < MSFeatureMap

  
  properties (SetAccess = protected)
    index      % Vector of indices used for picking output values
  end
  
  methods
    function obj = MSFeatureSelectionMap (numFeaturesIn, index, creator)
      
      narginchk(2,3);
      if nargin < 3
        % If creator string is not specified, use this class name
        creator = 'MSFeatureSelectionMap';
      end
      
      % Check input arguments
      if ~(~isempty(index) && isnumeric(index) && isvector(index))
          error('index must be an integer vector');
      end
      % Call super class constructor where mzVector validity is checked
      obj@MSFeatureMap([], creator, numFeaturesIn);
      % Check whether specified index values are unique and in proper range
      u = unique(index);
      if ~(u(1) >= 1 && u(end) <= numFeaturesIn && length(u) == length(index))
        error('index must contain unique index values in [1..%d]',numFeaturesIn);
      end
      obj.index = floor(index);
    end
  end
    
  methods (Access=protected)
    function featureData = map_impl (obj, msData, itemMask, numFeatures, ~)
      dataSelection = msData.data( :, obj.index(1:numFeatures) );
      dataSelection( ~itemMask, : ) = 0;
      featureData = MSFeatureData( dataSelection, class(obj), msData, itemMask );
      if isa( msData, 'MSMaldiData' )
        featureData.featureInfo.mzValues = msData.mzVector(obj.index(1:numFeatures));
      end
    end
    
    function numFeatures = getNumFeatures(obj)
      % Get number of output features  
      numFeatures = length(obj.index);
    end    
  end
  
end
  
  
  