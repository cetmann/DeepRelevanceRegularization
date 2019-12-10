classdef MSNeuralNetFeatureMap<MSFeatureMap
    % Generates feature data using the whole mz-range in msData. The
    % selected spectra remain unmodified in the corresponding feature data
    % rows.
    % Methods
    %  map_impl: generates feature data
    methods
        function obj=MSNeuralNetFeatureMap(mzVector, creator)
            obj@MSFeatureMap(mzVector,creator);
        end
    end
    methods(Access=protected)
        
        function featureData=map_impl(obj, msData, itemMask, numFeatures, ~)
            % generates feature data consisting on unmodified spectra
            % specified by itemMask
            featureData = MSFeatureData(msData.data);
        end
        
        function numFeatures = getNumFeatures (~)
            numFeatures = inf;
        end
    end
end