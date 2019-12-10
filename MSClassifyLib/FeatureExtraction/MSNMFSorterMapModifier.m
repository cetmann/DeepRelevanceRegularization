classdef MSNMFSorterMapModifier < MSMapModifier
    % Class for the application of basis sorting to a basis-based feature
    % map. The basis are sorted according to the sorting criterium 
    %
    % Properties
    %   sorting: sorting method (string)
    %
    % Methods
    %   MSNMFSorterMapModifier: Constructor
    %   apply: returns the modified feature map with sorted basis
    %   isValidSorting: returns whether the input sorting is valid (logical)
    %
    % MSNMFSorterMapModifier uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    
    properties (SetAccess = immutable)
        sorting; %sorting method (string)
    end
    methods
        function obj = MSNMFSorterMapModifier(sorting)
            %Constructor
            %obj=MSNMFSorterMapModifier-the default 'l1l2' is applied
            %obj=MSNMFSorterMapModifier(sorting)
            %the input string sorting is applied
            narginchk(0,1)
            if nargin<1
                obj.sorting = 'l1l2';
            else
                if ~MSNMFSorterMapModifier.isValidSorting(sorting)
                    error('unknown specified sorting')
                end
                obj.sorting=sorting;
            end
            obj.validMaps = {'MSBasisMap'};            
        end
        
        function featureMap=apply(obj, nmfMap, msData, msLabels)
            %obj.apply(nmfMap, msData, msLabels)
            %INPUT
            %nmfMap: MSFeatureMap object
            %msData: MSMaldiData object, necessary for roc-based sorting
            %msLabels: MSLabelData object, necessary for roc-based sorting
            %OUTPUT
            %featureMap: new modified (sorted) MSFeatureData
            %If the sorting is 'initial', 'l1l2' or 'gini' the msData and
            %msLabels arguments are ignored
            
            featureMap=MSBasisMap(nmfMap.decompositionObj.copy,nmfMap.mzVector,nmfMap.creator);
            featureMap.featureMapInfo=nmfMap.featureMapInfo;
            featureMap.switchProjectionType(nmfMap.projectionType);
            if any(strcmp(obj.sorting,{'initial','l1l2','gini'}))
                featureMap.decompositionObj.sortBasis(obj.sorting);
            else
                featureMap.decompositionObj.sortBasis(obj.sorting, msData, msLabels);
            end
        end
    end
    methods(Static)
        function bool = isValidSorting(sorting)
            %bool=MSNMFSorterMapModifier.isValidSorting(sorting)
            %Returns whether the input string sorting is a valid sorting
            bool=any(strcmp(sorting,{'initial','l1l2','gini','rocMax',...
                'rocMaxLabelBalance'}));  
        end
    end
   

end