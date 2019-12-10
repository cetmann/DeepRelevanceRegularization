classdef MSSupervisedFeatureSelection < MSFeatureExtraction
  % Feature extraction based on supervised mz-values selection
  %

  properties (SetAccess = protected)
    allowedDistanceMeasures
    allowedMulti2Binary
    allowedBinarySetAggregation
  end
  properties
    distanceMeasure
    multi2Binary
    binarySetAggregation
  end
  
  methods
    function obj = MSSupervisedFeatureSelection(maxNumFeatures, distanceMeasure, multi2Binary, binarySetAggregation )
      narginchk(0,4);
      
      % Initialize properties
      if nargin < 4
        binarySetAggregation = 'mean';
      end
      if nargin < 3
        multi2Binary = 'oneVsOne';
      end
      if nargin < 2
        distanceMeasure = 'bhattacharrya';
      end        
      if nargin < 1
        maxNumFeatures = inf;
      end
      obj@MSFeatureExtraction(maxNumFeatures);
      
      obj.allowedDistanceMeasures = {'l1','bhattacharrya','corr','feature_corr'};
      obj.allowedMulti2Binary = {'oneVsOne','oneVsRest'};
      obj.allowedBinarySetAggregation = {'min','max','median','mean'}; 

      obj.distanceMeasure = distanceMeasure;
      obj.multi2Binary = multi2Binary;
      obj.binarySetAggregation = binarySetAggregation;
    end
  end
  
  methods (Access=protected)
    function map = createMap_impl (obj, msData, labels)
      classes = unique(labels(labels~=0));
       nClasses = length(classes);
      if nClasses < 2
          error('Not enough labels specified to create the feature map')
      end
      
      indexSets = cell( 2, nClasses*(nClasses+1)/2 );
      k = 0;
      switch obj.multi2Binary
        case 'oneVsOne'
          for i=1:nClasses
            for j=i:nClasses
              l1 = labels == classes(i);
              l2 = labels == classes(j);
              k = k+1;
              indexSets{ 1, k } = l1;
              indexSets{ 2, k } = l2;
            end
          end 
        case 'oneVsRest'
          for i=1:nClasses
            l1 = labels == classes(i);
            l2 = ~l1 & labels ~= 0;
            k = k+1;
            indexSets{ 1, k } = l1;
            indexSets{ 2, k } = l2;            
          end
      end
      
      switch obj.binarySetAggregation
        case 'min'
          agg = @(x) min(x,[],2);
        case 'max'
          agg = @(x) max(x,[],2);
        case 'mean'
          agg = @(x) mean(x,2);
        case 'median'
          agg = @(x) median(x,2);
      end
      
      if obj.distanceMeasure == 'feature_corr'
            indexWeights = zeros(msData.dataLength,k);
            for i = 1:k
                l1 = indexSets{1,i};
                l2 = indexSets{2,i};

                m1 = msData.data( l1, : );
                m2 = msData.data( l2, : );
                idxs = find(abs(corrcoef(m1)-corrcoef(m2))>0.2);
                [I,~] = ind2sub([msData.dataLength msData.dataLength],idxs);
                indexWeights(:,i) = histc(I,(1:msData.dataLength))/msData.dataLength;
            end
            indexWeights = agg(indexWeights);
      else
          switch obj.distanceMeasure
            case 'l1'
              distance = @MSSupervisedFeatureSelection.dist_l1;
            case 'bhattacharrya'
              distance = @MSSupervisedFeatureSelection.dist_bhattacharrya;
            case 'corr'
              distance = @MSSupervisedFeatureSelection.dist_corr;
          end

          indexWeights = zeros(msData.dataLength,1);
          for j = 1:msData.dataLength
            dists = zeros(k,1);
            for i = 1:k
                l1 = indexSets{1,i};
                l2 = indexSets{2,i};

                mTotal = msData.data( l1|l2, j );
                m1 = min(mTotal);
                m2 = max(mTotal);

                N = ceil(min(sum(l1),sum(l2))/20);
                if N == 0
                  error('Empty class');
                end
                if N < 5
                  warning('One class has very few items');
                end
                x = linspace( m1, m2, N );

                h1 = hist( msData.data(l1,j), x);
                h2 = hist( msData.data(l2,j), x);

                dists(i) = distance( h1,h2 );
            end
            indexWeights(j) = agg(dists);
          end
      end

      [indexWeights,indexList] = sort(indexWeights,'descend');
      
      numFeatures_data = min( obj.maxNumFeatures, msData.dataLength );
      map = MSFeatureSelectionMap( msData.dataLength, ...
        indexList(1:numFeatures_data), class(obj) );                     
      map.featureMapInfo.indexWeights = indexWeights(1:numFeatures_data);
      map.featureMapInfo.rankedList = indexList;
      map.featureMapInfo.weights = indexWeights;
      
    end
  end
  
  methods(Static)
    function d = dist_l1(x,y)
      d = sum(abs( MSSupervisedFeatureSelection.norml1(x)-...
        MSSupervisedFeatureSelection.norml1(y) ));      
    end
    function d = dist_bhattacharrya(x,y)
      bc = sum( sqrt( MSSupervisedFeatureSelection.norml1(x).*...
        MSSupervisedFeatureSelection.norml1(y) ) );
      bc = min( bc, 1 );
      
      d = sqrt( 1 - bc );
    end
    function d = dist_corr(x,y)
      d = (1-corr(x(:),y(:)))/2;
    end
    
    
    function x = norml1(x)
      x = x/sum(x);
      assert( all( x>=0 ), 'negative values' );
    end
  end
  
  
  
end

