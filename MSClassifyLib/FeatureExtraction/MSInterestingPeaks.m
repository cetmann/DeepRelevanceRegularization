classdef MSInterestingPeaks < MSFeatureExtraction
  % Feature extraction based on unsupervised m/z-value selection
  
  methods
    function obj = MSInterestingPeaks (maxNumFeatures)
      if nargin < 1
        maxNumFeatures = inf;
      end
      % Initialize properties
      obj@MSFeatureExtraction(maxNumFeatures);
    end
  end
  methods (Access=protected)
    function map = createMap_impl (obj, msData, ~)
      spcMean = msData.meanData();
      spcMed = median(msData.data,1);
      spcMAD = median(abs(msData.data-repmat(spcMed,msData.numItems,1)),1);
      [~,idx1] = sort(spcMean,'descend');
      [~,idx2] = sort(spcMAD,'descend');
      idx12 = [idx1;idx2];
      [~,idx12uIdx] = unique(idx12);
      idxF = idx12( sort(idx12uIdx,'ascend') );
      map = MSMzValuesMap(msData.mzVector, ...
        idxF(1:obj.maxNumFeatures), class(obj) );
    end
  end
end





