classdef MSHybridFeatureExtraction < MSFeatureExtraction

    properties( SetAccess = immutable )
       mainFE
       addFE 
    end

    
  methods
    function obj = MSHybridFeatureExtraction ( CFE )
        narginchk(1,1)
        obj@MSFeatureExtraction( CFE{1}.maxNumFeatures );
        obj.mainFE = CFE{1};
        obj.addFE = CFE(2:end);
        obj.supportsDataPartition = true;
        
    end
  end
    
  methods ( Access=protected )
    function featureMap = createMap_impl (obj, msData, labels)    
        featureMap = obj.call_createMap(obj.mainFE, msData, labels );
        numAddD = length(msData.additionalDataVersions);
        numAddFE = length(obj.addFE);
        if numAddD ~= numAddFE
            warning('Number of add. data versions and add. FE does not match');
        end
        for i = 1:min(numAddD,numAddFE)
            featureMap.additionalFeatureMaps{i} = obj.call_createMap( ...
                obj.addFE{i}, msData.additionalDataVersions{i}, labels );
        end
    end
  end
    
    
    methods( Static ) 
        function featureMap = call_createMap ( FE, msData, labels )
           if FE.supportsDataPartition || ~isa( labels, 'MSDataPartition' )
               featureMap = FE.createMap( msData, labels );
           else
               featureMap = FE.createMap( msData, labels.classes );
           end            
        end
    end
   
end