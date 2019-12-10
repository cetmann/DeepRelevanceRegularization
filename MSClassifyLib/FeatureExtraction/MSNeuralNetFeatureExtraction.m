classdef MSNeuralNetFeatureExtraction<MSFeatureExtraction
    % Dummy feature extraction class used to create a dummy feature map,
    % which leaves the received maldi data unmodified as feature data
    % Methods
    %   createMap_impl: Creates the dummy map
    methods (Access=protected)
        
        function map = createMap_impl(obj,msData,~)
            % creates dummy map
            map=MSNeuralNetFeatureMap(msData.mzVector, class(obj));
        end
    end
end