classdef MSMixedMapModifier < MSMapModifier
    % Class implementing a map modifier for the class MSMixedMapModifier.
    % Initialized with a list of map modifications it performes those
    % modifications for each map constituting the mixed map (if possible)
    % in the same order that the list of modifications is given.
    %
    % properties:
    % listMM: cell array of map modifications
    % nMM (dependent): number of map modifications.
    %
    % methods
    % MSMixedMapModifier: constructor
    % apply: Applies secuence of map modifiers to mixed feature map
    
    properties (SetAccess = immutable)
        listMM; % list of map modifiers
    end
    properties(Dependent)
        nMM; % number of map modifiers
    end
    
    methods
        function obj = MSMixedMapModifier(listMapMod)
            % obj = MSMixedMapModifier (Constructor)
            % INPUT
            % listMapMod: List of map modifications or a map modifier
            % object
            listMapMod = assert(listMapMod);
            obj.listMM = listMapMod;
            % obj.validMaps = {'MSMixedFeatureMap'};
        end
        function numMM = get.nMM(obj)
            numMM = length(obj.listMM);
        end

        function featureMap=apply(obj, inputFeatureMap, msData, msLabels)
            if isa(inputFeatureMap, 'MSMixedFeatureMap')
                listFM = inputFeatureMap.listFM;
            else
                listFM = {inputFeatureMap};
            end
            % apply each map modification to each possible map in list of
            % maps
            for iMM = 1:obj.nMM
                mapMod = obj.listMM{iMM};
                for iFM = 1:length(listFM)
                    map = listFM{iFM};
                    if mapMod.applicable(map)
                        listFM{iFM} = mapMod.apply(map, msData, msLabels);
                    end
                end
            end
            if length(listFM )> 1
                featureMap = MSMixedFeatureMap (inputFeatureMap.mzVector, listFM,...
                                            inputFeatureMap.creator);
            else
                featureMap = listFM{1};
            end
        end
    end
    
end
function listMM = assert(listMapMod)
    if (isa(listMapMod,'MSMapModifier'))
        listMM = {listMapMod};
        return
    end
    if ~iscell(listMapMod)
        error('The input must be either a MapModifier or a list of them');
    end
    for i=1:length(listMapMod)
        if ~isa(listMapMod{i},'MSMapModifier')
            error('The input must be either a MapModifier or a list of them');
        end
    end
    listMM = listMapMod;
end