classdef MSHardThrMapModifier < MSMapModifier
    % Class for the application of hard thresholding to feature map basis. 
    % A threshold (property of the class) is applied to each basis vector
    % setting to zero those values bellow that value.
    % The threshold value is specified as a percentile value, e.g. thr ==
    % 99.5 means only the 0.5% highest peaks will keep their height while
    % the rest is set to 0
    %
    % Properties
    %   threshold: threshold value
    %
    % Methods
    %   MSHardThrMapModifier: Constructor
    %   apply: returns the modified feature map
    %
    % MSHardThrMapModifier uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    
    properties (SetAccess=immutable)
        threshold; % threshold value (positive real)
    end
    methods
        function obj=MSHardThrMapModifier(threshold)
            %Constructor
            %obj=MSHardThrMapModifer(threshold)
            %INPUT
            %threshold: threshold value (positive real)
            if ~(isscalar(threshold)&&isnumeric(threshold)&&threshold>0)
                    error('threshold must be a positive number');
            end
            obj.threshold=threshold;
            obj.validMaps = {'MSBasisMap'};
        end
        function nmfThrMap=apply(obj, nmfMap,~,~)
            %obj.apply(nmfMap,~,~)
            %Returns a new MSFeatureMap nmfThrMap by applying the threshold
            %to the input MSFeatureMap nmfMap
           nmfThrMap = MSBasisMap(nmfMap.decompositionObj.copy,nmfMap.mzVector, nmfMap.creator);
           nmfThrMap.decompositionObj.unsortedBasis=MSHardThr(nmfMap.decompositionObj.unsortedBasis,obj.threshold);
           nmfThrMap.featureMapInfo=nmfMap.featureMapInfo;
           nmfThrMap.switchProjectionType(nmfMap.projectionType);
           nmfThrMap.decompositionObj.sortTypeCurr = nmfMap.decompositionObj.sortTypeCurr;
        end 
    end
end