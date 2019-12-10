classdef MSMixedFeatureData < MSFeatureData
    % Defines a type of feature data containing features generated from
    % different maps.
    %
    % Properties
    %   numFeaturesPerType: array with number of features per map type
    %   types: cell array with class names of the feature maps generating
    %   the different type of features
    properties (SetAccess = private)
        numFeaturesPerType;
        types;
    end
    methods
        function obj = MSMixedFeatureData(numFeaturesPerType, types, data, creator, sourceData, validItems)
            narginchk(3,6)
            if nargin < 4
                creator = '';
            end
            if nargin < 5
                sourceData = [];
            end
            if nargin < 6
                validItems = [];
            end
            obj@MSFeatureData(data, creator, sourceData, validItems)
            % check first 2 input arguments
            if ~(isvector(numFeaturesPerType)&&isnumeric(numFeaturesPerType)...
                &&all(numFeaturesPerType>0)&&all(floor(numFeaturesPerType)==numFeaturesPerType))
                error('The first argument must be a vector of natural numbers')
            end
            if sum(numFeaturesPerType) ~= obj.dataLength
                error('The sum of numbers of features per type must be equal to the dataLength');
            end
            if ~iscellstr(types)||length(types)~=length(numFeaturesPerType)
                error('The number of elements in first and second input parameters must be equal');
            end
            obj.numFeaturesPerType = numFeaturesPerType;
            obj.types = types;
        end
        function fd = reduceFeatures(obj, numFeatures)
            % Returns a feature data object with a reduced number of
            % features. Features are selected proportionally to the number
            % of features for each type of map in the feature definition.
            % fd = obj.reduceFeatures(numFeatures)
            % INPUT
            %   numFeatures: number of resulting features
            % OUTPUT
            %   fd: feature data with reduced number of features
            if numFeatures > obj.dataLength
               % warning('Nothing to reduce. The actual number of features is not larger than input')
            return
            end
            % get number of output features per feature type
            nFeaturesList = MSMixedFeatureData.nFeaturesPerTypeOut(numFeatures,obj.numFeaturesPerType);
            data = zeros(obj.numItems, numFeatures);
            countTotal = 0;
            count = 0;
            % build data according to number of features per type
            for i = 1:length(obj.numFeaturesPerType)
                data(:,count + 1:count + nFeaturesList(i))=obj.data(:,countTotal + 1:countTotal+nFeaturesList(i));
                count = count + nFeaturesList(i);
                countTotal = countTotal + obj.numFeaturesPerType(i);
            end
            mask = nFeaturesList > 0;
            if sum(mask) == 1 % Generate a simple feature data if only one type of features present
                type = obj.types{mask};
                fd = MSFeatureData(data, type, obj, obj.validItems);
            else      % if not, generate a mixed feature map.
                fd = MSMixedFeatureData(nFeaturesList(mask), obj.types(mask),...
                                data, obj.creator, obj, obj.validItems);
            end
        end
    end
    methods(Static)
        function nFeaturesOut = nFeaturesPerTypeOut (numFeatures, numFeaturesPerType)
            % Returns an array with the number of features per type, given
            % that a certain number of total features must be selected.
            % This is done proportionally to the number of features per
            % type.
            % nFeaturesOut = MSMixedFeatureData.nFeaturesPerTypeOut (numFeatures, numFeaturesPerType)
            % INPUT:
            %   numFeatures: number of features to be selected
            %   numFeaturesPerType: number of features per type available
            % OUTPUT:
            %   nFeaturesOut: number of features per type selected.
            total = sum(numFeaturesPerType);
            nFeaturesOut = floor(numFeaturesPerType/total*numFeatures);
            toplace = numFeatures - sum(nFeaturesOut);
            [~,index] = sort(numFeaturesPerType,'descend');
            nFeaturesOut(index(1:toplace)) = nFeaturesOut(index(1:toplace))+1;
        end
    end
end