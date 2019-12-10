classdef MSBaseLineCorrMapModifier<MSMapModifier
    % Class for the baseline correction of basis vectors of an
    % MSLinearBasisMap map.
    %
    % Methods
    %   MSBaseLineCorrMapModifier: Constructor
    %   apply: Returns a new map from the input map with the baseline
    %          correction
    %
    % MSBaseLineCorrMapModifier uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    
    methods
        function obj=MSBaseLineCorrMapModifier()
            obj.validMaps = {'MSBasisMap'};
        end
        function nmfBaseMap=apply(~, nmfMap,msData,~)
            %obj.apply(nmfMap,msData,~)
            %INPUT
            %nmfMap: MSFeatureMap to be modified
            %msData: MSMaldiData object
            %Returns a new map from the input map with the baseline
            %correction. It uses the msbackadj matlab function
            nmfBasis = abs(msbackadj(msData.mzVector', nmfMap.decompositionObj.basis'))';
            decomposition=nmfMap.decompositionObj.copy;
            decomposition.basis=nmfBasis;
            nmfBaseMap = MSBasisMap(decomposition, nmfMap.mzVector, nmfMap.creator);
            nmfBaseMap.featureMapInfo = nmfMap.featureMapInfo;
            nmfBaseMap.featureMapInfo.numNonZeros = sum(logical(...
                                    nmfBaseMap.decompositionObj.basis),2);
            nmfBaseMap.switchProjectionType(nmfMap.projectionType);
        end 
    end
end