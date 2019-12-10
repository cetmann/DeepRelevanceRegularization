classdef MSMixedFeatureMap < MSFeatureMap
    % Defines a mixed feature map class where different features (e.g., ROC
    % and NMF are contained consecutively in the same map)
    % Properties
    %   listFM: list of feature map objects
    %   numFeaturesTotal: max number of combined features
    %   numFeaturesPerFM: number of features per map type.
    %   nFM: number of feature maps
    %
    % Methods
    %   MSMixedFeatureMap: Constructor
    
    properties (SetAccess = immutable)
        listFM;
        numFeaturesTotal;
        numFeaturesPerFM;
    end
    properties (Dependent)
        nFM;        
    end
    methods
        function obj = MSMixedFeatureMap(mzVector, listFM, creator)
            narginchk(2,3)
            if nargin < 3
                creator = '';
            end
            obj@MSFeatureMap(mzVector, creator);
            assert(listFM);
            obj.listFM = listFM;
            obj.setNumFeatures;
            [obj.numFeaturesTotal, obj.numFeaturesPerFM] = obj.setNumFeatures;
        end
        function numFM = get.nFM(obj)
            numFM = length(obj.listFM);
        end

    end
    methods (Access = protected)
        function featureData = map_impl(obj, msData, itemMask, numFeatures, printFlag)
            data = [];
            dataType = cell(1,obj.nFM);
                     
            for i=1:obj.nFM
                % get feature data per feature map
                fd = obj.listFM{i}.map(msData, itemMask, obj.listFM{i}.numFeatures, printFlag);
                % aggregate feature data to feature data already generated
                % from other maps
                data =[data fd.data;];                
                % store feature map type
                dataType{i} = class(obj.listFM{i});
            end
 
            featureData = MSMixedFeatureData(obj.numFeaturesPerFM, dataType, data, class(obj), msData, itemMask);
            featureData = featureData.reduceFeatures(numFeatures);
        end
        
        function numFeatures = getNumFeatures(obj)
            numFeatures = 0;
            for i=1:obj.nFM
                numFeatures = numFeatures + obj.listFM{i}.numFeatures;
            end
        end
    end
    
    methods (Access = private)
        function [numFeaturesTotal, numFeaturesPerFM] = setNumFeatures(obj)
            % used to set properties obj.numFeaturesTotal and
            % obj.numFeaturesPerFM
            numFeaturesTotal = obj.getNumFeatures;
            numFeaturesPerFM = zeros(1,obj.nFM);
            for i = 1:obj.nFM
                nf = obj.listFM{i}.numFeatures;
                numFeaturesPerFM(i) = nf;                
            end
        end
        
        function numFeaturesArray = distributeNumFeatures(obj, numFeatures)
            [~, index] = sort(obj.numFeaturesPerFM, 'ascend');
            if numFeatures > obj.numFeaturesTotal
                numFeaturesArray = obj.numFeaturesPerFM;
                return;
            end
            numFeaturesArray = zeros(1, obj.nFM);
            % - Distribute num of fd per type proportionally to initial
            %   distribution
            features = floor((numFeatures/obj.numFeaturesTotal)*(obj.numFeaturesPerFM));                
            numFeaturesArray = numFeaturesArray + features;
            numFeatures = numFeatures - sum(features);
            % - Asign rest of features according to priority order
            numFeaturesArray(index(1:numFeatures)) = numFeaturesArray(index(1:numFeatures))+1;          
            
        end
    end

end
function assert (listFM)
     errormsg = ['The input argument must be an MSFeatureExtraction' ...
                'object or a list of feature extractions'];

    if ~(iscell(listFM) && numel(listFM)>1)
        error(errormsg)
    else
        for i=1:length(listFM)
            if ~isa(listFM{i},'MSFeatureMap')
                error(errormsg)
            end
        end
    end
end