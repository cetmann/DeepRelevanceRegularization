classdef MSSelectMapsMapModifier < MSMapModifier
    % class defining "select maps" map modification which consists on
    % selecting a subgroup of map types from a mixed map. This way one
    % could generate through the apply method an MSMzValuesMap (ROC based)
    % or an MSBasisMap (NMF based) from a mixed map with both maps.
    %
    % Properties
    %   selection: cell string with names of map type classes that should
    %   be selected.
    properties(SetAccess = immutable)
        selection
    end
    methods
        function obj = MSSelectMapsMapModifier(selection)
            selection = assert(selection);
            obj.selection = selection;
            obj.validMaps = {'MSMixedFeatureMap'};
        end
        function featureMap=apply(obj, inputFeatureMap, ~, ~)
            % find mask indicating which elements of the list of maps of
            % the input mixed map will be selected
            indexMask = obj.indexSelection(inputFeatureMap);  
            % Get corresponding maps
            listMaps = inputFeatureMap.listFM(indexMask);
            switch length(listMaps)
                case 0 % No map could be selected-> create dummy map
                    featureMap = MSDummyFeatureMap(obj.mzVector, obj.creator);
                    warning('The resulting modified map is a dummy map')
                case 1 % Only one map selected-> Results in a simple feature map object
                    featureMap = listMaps{1}.copy;
                otherwise % More than one map-> Results in a mixed feature map
                    featureMap = MSMixedFeatureMap(inputFeatureMap.mzVector,...
                                      listMaps, inputFeatureMap.creator);
            end
        end
    end
    methods(Access = private)
        function indexMask = indexSelection(obj, featureMap)
            sel = obj.selection;
            lFM = featureMap.listFM;
            indexMask = false(1, length(lFM));
            for iFM =1:length(lFM)
                for iSel = 1:length(sel)
                    if isa(lFM{iFM},sel{iSel})
                        indexMask(iFM) = true;
                        break;
                    end
                end
            end
        end
    end
    
end
function selectionOut = assert(selection)
    if ischar(selection)
        selection = {selection};
    elseif ~iscellstr(selection)
        error('the input selection must be a cell array of strings containing name of class type')
    end
    % check if strings correspond to MSFeatureMap class or subclasses.
    for i = 1:length(selection)
        try
            name = meta.class.fromName(selection{i}).SuperclassList.Name;
        catch
            error('The input contains elements which are not a valid FeatureMap name')
        end
        if ~strcmp(name, 'MSFeatureMap')
            error('each string of the input cell must be a valid MSFeatureMap subclass name')
        end
    end
    selectionOut = selection;
end