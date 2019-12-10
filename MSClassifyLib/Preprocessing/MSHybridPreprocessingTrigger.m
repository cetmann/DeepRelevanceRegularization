classdef MSHybridPreprocessingTrigger < MSPreprocessingTrigger

    properties( SetAccess = immutable )
        mainCPP;
        addCPP;
    end

    
  methods
    function obj = MSHybridPreprocessingTrigger ( CPP )
        narginchk(1,1)
        obj.mainCPP = CPP{1};
        obj.addCPP = CPP(2:end);
    end
    
    function apply(obj, msData)
        numAddD = length(msData.additionalDataVersions);
        numAddCPP = length(obj.addCPP);
        for i = 1:numAddCPP
            if i > numAddD 
                msData.additionalDataVersions{i} = msData.copy;
            end
            obj.addCPP{i}.apply(msData.additionalDataVersions{i});
        end
        obj.mainCPP.apply(msData);
    end
    
  end
   
end