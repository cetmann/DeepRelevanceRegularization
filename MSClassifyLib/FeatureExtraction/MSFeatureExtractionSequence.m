classdef MSFeatureExtractionSequence < MSFeatureExtraction
    % Defines a feature extraction method concatenating different types of
    % feature extractions in the sense of subsequent mappings from left to
    % right
    % (e.g. ROCPicking followed by NMF on the list of ROC values)
    % 
    %   Properties
    %   listFE: List of feature extractions (cell)
    %   nFE (Dependent): Number of feature extractions in definition
    %   nFOut (Dependent): Number of features extracted in total 
    %
    %   Methods
    %   MSFeatureExtractionSequence: Constructor
    properties (SetAccess = immutable)
        listFE = [];
        % maxNumFeatures in superclass
    end
    properties (Dependent)
        nFE;
        nFIntermediate;
    end
    methods
        function obj = MSFeatureExtractionSequence(listFE)
            % obj = MSMixedFeatureExtraction: Constructor
            % Defines a feature extraction method mixing different types of
            % features
            % INPUT
            %   listFE: list of feature extraction objects
            
            nFE = length(listFE);
            nFOut = listFE{nFE}.maxNumFeatures;
            obj@MSFeatureExtraction( nFOut );
            obj.listFE = listFE;
            MSFeatureExtractionSequence.assert(obj);
        end
        function numFE = get.nFE(obj)
            numFE = length(obj.listFE);
        end
        function nFIntermediate = get.nFIntermediate(obj)
            nnFE = obj.nFE;
            nFIntermediate = zeros(1,nnFE);
            for k = 1:nnFE
                nFIntermediate(k) = obj.listFE{k}.maxNumFeatures;
            end
        end
    end
    methods (Access=protected)
        function featureMap = createMap_impl(obj,msData,labels)
  
            listFM = cell(1, obj.nFE);
            % Creates individual feature maps and generates feature map
            % sequence
            
            fData = msData;
            for k = 1:obj.nFE
                listFM{k} = obj.listFE{k}.createMap(fData,labels);
                if k < obj.nFE
                    fData = listFM{k}.map(fData);
                    dummyMzVector = 1:fData.dataLength;
                    fData = MSMaldiData( fData.data, dummyMzVector );
                end
            end
            featureMap = MSFeatureMapSequence(msData.mzVector, listFM, class(obj));
               
        end
    end
    methods (Static)
        function assert(obj)
            errormsg = ['The input argument must be a list of feature '...
                        'extractions with at least 2 elements'];
            if ~(iscell(obj.listFE) && numel(obj.listFE)>1)
                error(errormsg)
            else
                for k=1:length(obj.listFE)
                    if ~isa(obj.listFE{k},'MSFeatureExtraction')
                        error(errormsg)
                    end
                end
            end
            warningmsg = ['Subsequent Feature extractions .maxNumFeatures'...
              ' properties are not descending'];
            nFIntermediate = obj.nFIntermediate;
            nFIDiff = diff( nFIntermediate );
            if any( nFIDiff > 0 )
                warning(warningmsg);
            end
        end
    end
end
