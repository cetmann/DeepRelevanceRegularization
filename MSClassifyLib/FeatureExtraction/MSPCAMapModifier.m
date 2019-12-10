classdef MSPCAMapModifier < MSMapModifier
  % Class for the PCA map modification. It combines the input MSFeatureMap in
  % apply with the MSPCAProjection property to generate a new map
  % map.
  %
  % Properties
  %  pCAProjection: MSPCAProjection object
  %
  % Methods
  %  MSPCAMapModifier: Constructor
  %  apply: Returns an MSMapPlusPCAMap combining a PCA projection and a feature
  %         map
  %
  % MSPCAMapModifier uses the handle semantic, i.e. when 
  % assigning an object of this class to a variable, only a reference to
  % the original object is copied. Use the copy method to create a deep copy.
    properties (SetAccess = immutable)
        pCAProjection %MSPCAProjection object
    end
    methods
        function obj = MSPCAMapModifier(numFeatures)
            %Constructor
            %obj=MSPCAMapModifier(numFeatures)
            %INPUT
            %numFeatures: number of features for the PCA projection
            %           property
            
            %input validation
            narginchk(0,1)
            if nargin<1
                obj.msPCAProjection=MSPCAProjection();
            else
                if ~isnumeric(numFeatures)||~isscalar(numFeatures)||~(numFeatures>=1)
                    error('numFeatures must be inf or a positive integer')
                end
                obj.pCAProjection=MSPCAProjection(floor(numFeatures));
            end            
            obj.validMaps = [];            
        end
        
        function featureMap=apply(obj, featureMap, ~, ~)
            %featureMap=obj.apply(mzMap,~,~)
            %Returns a feature map combining the PCA
            %projection property with the input map
            %INPUT
            % featureMap: MSFeatureMap object
            %OUTPUT
            % featureMap: MSMapPlusPCAMap object
            
            featureMap = MSMapPlusPCAMap(featureMap.mzVector, obj.pCAProjection, featureMap); 
            
        end
    end
end