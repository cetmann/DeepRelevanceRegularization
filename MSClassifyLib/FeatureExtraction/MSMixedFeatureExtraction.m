classdef MSMixedFeatureExtraction < MSFeatureExtraction
    % Defines a feature extraction method mixing different types of features
    % (e.g. NMF and ROC features)
    % 
    %   Properties
    %   listFE: List of feature extractions
    %   nFE (Dependent): number of feature extractions in definition
    %
    %   Methods
    %   MSMixedFeatureExtraction: Constructor
    properties (SetAccess = immutable)
        listFE = [];
    end
    properties (Dependent)
        nFE;
    end
    methods
        function obj = MSMixedFeatureExtraction(listFE)
            % obj = MSMixedFeatureExtraction: Constructor
            % Defines a feature extraction method mixing different types of
            % features
            % INPUT
            %   listFE: list of feature extraction objects
            
            obj.assert(listFE);          
            obj.listFE = listFE;
        end
        function numFE = get.nFE(obj)
            numFE = length(obj.listFE);
        end
    end
    methods (Access=protected)
        function featureMap = createMap_impl(obj,msData,labels)
            listFM = cell(1, obj.nFE);
            % creates individual feature maps and generates mixed feature map 
            for i = 1:obj.nFE
                listFM{i} = obj.listFE{i}.createMap(msData,labels);
            end
            featureMap = MSMixedFeatureMap(msData.mzVector, listFM, class(obj));
        end
    end
    methods (Static)
        function assert(listFE)
            errormsg = ['The input argument must be a list of feature '...
                        'extractions with at least 2 elements'];
            if ~(iscell(listFE) && numel(listFE)>1)
                error(errormsg)
            else
                for i=1:length(listFE)
                    if ~isa(listFE{i},'MSFeatureExtraction')
                        error(errormsg)
                    end
                end
            end
        end
    end
end
