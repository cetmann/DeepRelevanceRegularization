classdef MSDummyMapModifier < MSMapModifier
    % Dummy class for map modification which does not modify the feature
    % map. This class is necessary to build a generic classification
    % pipeline and is used when no modifications of the original feature
    % map are required
    % Methods
    %   apply: returns the same feature map
    %
    % MSDummyMapModifier uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    methods
        function featureMap=apply(~, inputFeatureMap, ~, ~)
            %obj.apply(inputFeatureMap)
            %returns exactly the same inputFeatureMap (the same reference)
            featureMap=inputFeatureMap.copy;
        end
    end
end