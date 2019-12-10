classdef MSFeatureMapSequence < MSFeatureMap
    % Defines a feature map sequence class where different features maps,
    % e.g. ROCPicking and NMF can be executed consecutively
    % Properties
    %   listFM: list of feature map objects
    %   numFeatures: output number of features after applying the last map
    %   nFM: number of feature maps
    %
    % Methods
    %   MSFeatureMapSequence: Constructor
    
    properties
        listFM;
    end
    properties (Dependent)
        nFM;
    end
    methods
        function obj = MSFeatureMapSequence(mzVector, listFM, creator)
            narginchk(2,3)
            if nargin < 3
                creator = '';
            end
            obj@MSFeatureMap(mzVector, creator);
            assert(listFM);
            obj.listFM = listFM;
        end
        function numFM = get.nFM(obj)
            numFM = length(obj.listFM);
        end

    end
    methods (Access = protected)

        function fData = map_impl(obj, msData, itemMask, numFeatures, printFlag)

              fData = msData;
              for k = 1:obj.nFM
                  if k < obj.nFM
                      fData = obj.listFM{k}.map(fData, itemMask, [], printFlag);
                      dummyMzVector = 1:fData.dataLength;
                      fData = MSMaldiData( fData.data, dummyMzVector );
                  else
                      fData = obj.listFM{k}.map(fData, itemMask, numFeatures, printFlag);
                  end
              end


        end
        
        function numFeatures = getNumFeatures(obj)
            numFeatures = obj.listFM{obj.nFM}.numFeatures;
        end  
      
      
    end
    
end
function assert (listFM)
     errormsg = 'The input argument must be a list of n>2 MSFeatureMaps';
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