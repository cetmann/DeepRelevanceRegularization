classdef MSPCAProjection < MSFeatureExtraction
  % Feature extraction based on principal component analysis
  %
  % This class implements a feature extraction based on a principal
  % component analysis (PCA) of the input data and a projection onto the
  % linear subspace defined by the principal components basis vectors.
  %
  %
  % Methods:
  %   MSPCAProjection: Constructor
  %   createMap_impl: Create MSLinearBasisMap from input data
  %
  % MSPCAProjection uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  methods
    function obj = MSPCAProjection (maxNumFeatures)
      % Constructor
      % obj = MSPCAProjection(numFeatures): Create PCA projection object
      %   numFeatures: Limit number of input features (optional)
      if nargin<1
          maxNumFeatures=inf;
      end
      obj@MSFeatureExtraction(maxNumFeatures);

    end
  end
  
  methods (Access=protected)
    
    function [map,d1,d2] = createMap_impl (obj, msData, labels)
      % Create linear basis map from PCA components of msData
      % map = obj.createMap(msData, itemMask):
      %   msData: Input dataset (MSData)
      %   labels: Optional mask vector selecting a subset of data
      %     items to be used for performing the PCA
      %
      %   The created map object is an MSLinearBasisMap, the featureMapInfo
      %   property contains a copy of the itemMask argument, the list of 
      %   eigenvalues and the computed explained variance for the PCA.
      
      itemMask=labels~=0;
      % Compute PCA
      N = min(msData.dataLength,obj.maxNumFeatures);
      [coeff,~,ev,~,explained] = pca(msData.data(itemMask,1:N), 'Centered',false);
      % Create feature map
      map = MSLinearBasisMap(msData.mzVector, coeff.', class(obj));
      % Store item mask, eigenvalues, explained variance in feature map info
      map.featureMapInfo.itemMask = itemMask;
      map.featureMapInfo.eigenvalues = ev;
      map.featureMapInfo.explainedVariance = explained;
      d1=[];
      d2=[];
    end
  end
  
end

