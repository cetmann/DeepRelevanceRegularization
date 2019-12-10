classdef MSChangeProjectionMapModifier<MSMapModifier
    properties (SetAccess=immutable)
        projectionType;
    end
    methods
        function obj = MSChangeProjectionMapModifier(projectionType)
            if ~ischar(projectionType)
                error('projectionType must be a string')
            end
            if ~strcmp(projectionType,MSBasisMap.listProjectionTypes)
                error('unknown projection type')
            end
            obj.projectionType=projectionType;
            obj.validMaps={'MSBasisMap'};
        end
        function featureMap=apply(obj, nmfMap, ~, ~)
            featureMap=nmfMap.copy;
            featureMap.switchProjectionType(obj.projectionType);            
        end        
    end
end