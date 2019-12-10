classdef MSMapModifier < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    % Base class for the map modifiers hierarchy. The functionality of this
    % class is to modify a given feature map. For example, when sorting the
    % basis of an NMF map, or removing the baseline, or applying PCA to a 
    % Roc map. This modification creates a new feature map and is performed
    % by the method <apply>
    %
    % Properties
    %   validMaps: Cell array of strings containing the names of the
    %         MSFeatureMap classes to which the modification can be applied
    %
    % Methods
    %   applicable: returns true if the modification can be applied  to the
    %           input map 
    %   apply (Abstract): returns a modified map
    %
    % MSMapModifier uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    properties
         validMaps=[]; %Cell array of strings containing the names of the
                        %SFeatureMap classes to which the modification can
                        %be applied
    end
    methods
        function bool = applicable(obj, featureMap)
         %obj.applicable(featureMap)
         %returns true if the modification can be applied  to the input map
         %Input:  MSFeatureMap featureMap
         %Output: logical bool
         
            if isempty(obj.validMaps)
                bool=isa(featureMap, 'MSFeatureMap');
            else
                %bool=any(strcmp(class(featureMap), obj.validMaps));
                
                % Necesary to do it like shown bellow if we admit that the 
                % validMap string represents a class up in the hierarchy. 
                % Desirable if in the feature the tree of MSFeatureData 
                % classes grows in depth
                
                for i=1:length(obj.validMaps)
                    if isa(featureMap,obj.validMaps{i})
                        bool=true;
                        return
                    end
                end
                bool=false;
            end
        end
    end
    methods (Abstract)
        % given a feature map returns a modified map
        featureMap=apply(obj, inputFeatureMap, msData, msLabels);
    end
end