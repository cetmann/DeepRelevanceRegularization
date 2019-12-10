classdef MSFeatureExtractionNMF < MSFeatureExtraction
  % Feature extraction based on non-negative matrix factorization with
  % the Kullback-Leibler discrepancy terms and in total 4 l1/l2 penalty 
  % terms
  %
  % This class implements a feature extraction based on a non-negative
  % matrix factorization (NMF) of the input data and a projection onto the
  % linear subspace defined by the NMF basis vectors.
  %
  %
  % Methods:
  %   MSFeatureExtractionNMFKL: Constructor
  %   createMap_impl: Create MSLinearBasisMap from input data
  %
  % MSNMFKLProjection uses the handle semantic, i.e. when assigning an 
  % object of this class to a variable, only a reference to the original 
  % object is copied. Use the copy method to create a deep copy.
  
  properties
      decomposerType='KL';%default (can be also orthoSpectra)
      mapParams=struct;
  end
  properties (Constant)
      listDecomposerTypes={'KL','orthoSpectra','frob','supScalProd',...
        'supScalProdOrthoSpectra','supLogReg','supADMM','targeted','triOrthoH'};
  end
  methods
    function obj = MSFeatureExtractionNMF (maxNumFeatures, decomposerType, varargin)
      % Constructor
      % obj = MSFeatureExtractionNMF(numFeatures): Create NMF projection object
      %   numFeatures: Define number of NMF features (default is set to 50)
      narginchk(0,inf)
      if nargin<1
          numFeatures = 50; % Default value
      else
          if~isscalar(maxNumFeatures)||~isnumeric(maxNumFeatures)||maxNumFeatures<=0
              error('maxNumber of features must be a positive integer')
          end
          numFeatures = floor(maxNumFeatures);
      end
      obj@MSFeatureExtraction(numFeatures);
      if nargin>=2 % get decomposerType
          if ~ischar(decomposerType)||~isvector(decomposerType)||...
              ~any(strcmp(decomposerType,obj.listDecomposerTypes))
              error('invalid decomposer type')
          end
          obj.decomposerType=decomposerType;
      end      
      %Extra map parameters, including the projection type;
      params=inputParser;
      params.KeepUnmatched = true;
      params.addParameter('checkT', 5, @isPosInteger); % override decomposer default
      params.addParameter('stopLim', 0.5*10^(-4), @isPosNumber); % override decomposer default
      params.addParameter('projectionType','correlation',...
          @(x) ischar(x)&& any(strcmp(MSBasisMap.listProjectionTypes,x)));
      params.addParameter('patterns',[],@validatePattern);

      %parse actual parameters
      params.parse(varargin{:});
%       for i=1:(length(varargin)-1)/2
%           obj.mapParams.(varargin{2*i})=varargin{2*i+1};
%       end
      %update properties    
      obj.mapParams=MSMergeAndTranslateDecStructs(params.Results, params.Unmatched);
      obj.supportLabelData = isSupervisedDecomposerType(obj.decomposerType);
    end
  end
  
  methods (Access=protected)    
    function map = createMap_impl (obj,msData,labels)
      % Create basis map out of msData via NMF
      % map = obj.createMap(msData, itemMask):
      %   msData: Input dataset (MSData)
      %   labels: non-supervised: Optional mask vector selecting a subset of data
      %     items to be used for performing the NMF,
      %           supervised: Required labelData for training
      %           (MSLabelData), also discarding items with label 0
      %      
      %   currently the default projection type (map type): 
      %           scalarProduct (can be changed by calling
      %             map.switchProjectionType)
      %   The created map object is an MSBasisMap
      
      % Check input arguments
      narginchk(2,3)
      supervised = isSupervisedDecomposerType(obj.decomposerType);
      % create parameters for decomposer 
      % (include modifications for new decomposer types in the method)
      varargin = obj.paramsForSpecificDecomposition(obj.mapParams, obj.decomposerType);
      %Decide which kind of decomposer will be used to generate the map
      switch obj.decomposerType
          case 'KL'% Choose MSDecomposerNMFKL
              decomposer = MSDecomposerNMFKL(obj.maxNumFeatures,varargin{:});
          case 'KLIsotope'% Choose MSDecomposerNMFKLIsotope
              decomposer = MSDecomposerNMFKLIsotope(obj.maxNumFeatures, msData.mzVector, varargin{:});
          case 'orthoSpectra' % Choose MSDecomposerNMFOrthoSpectra
              decomposer = MSDecomposerNMFOrthoSpectra(obj.maxNumFeatures,varargin{:});
          case 'frob' % Choose MSDecomposerNMFFrob
              decomposer = MSDecomposerNMFFrob(obj.maxNumFeatures,varargin{:});
          case 'frobIsotope' % Choose MSDecomposerNMFFrob
              decomposer = MSDecomposerNMFFrobIsotope(obj.maxNumFeatures, msData.mzVector, varargin{:});
          case 'supScalProd' % Choose MSDecomposerNMFSupScalProd
              decomposer = MSDecomposerNMFSupScalProd(obj.maxNumFeatures,varargin{:});
          case 'supScalProdOrthoSpectra' % Choose MSDecomposerNMFOrthoSupScalProd
              decomposer = MSDecomposerNMFOrthoSupScalProd(obj.maxNumFeatures,varargin{:});
          case 'supLogReg' % Choose MSDecomposerNMFSupervisedLogReg
              decomposer = MSDecomposerNMFSupervisedLogReg(obj.maxNumFeatures,varargin{:});
          case 'supADMM' % Choose MSDecomposerNMFSupADMM
              decomposer = MSDecomposerNMFSupADMM(obj.maxNumFeatures,varargin{:});
          case 'targeted'
              decomposer = MSDecomposerNMFTargeted(obj.maxNumFeatures, obj.mapParams.patterns,...
                                                   varargin{:});
          case 'triOrthoH'
              decomposer = MSDecomposerNMFTriOrthoH(obj.maxNumFeatures, varargin{:});          
      end
      % check labels
      if supervised
        % supervised decomposer type
        if nargin < 3 || ~isa(labels, 'MSLabelData')
          error('MSLabelData required as third argument for supervised decomposer type');
        end
        itemMask = labels.data~=0;
      else
        % non-supervised decomposer type
        if nargin >= 3 && ~isempty(labels)
          itemMask = labels~=0; 
        else
          itemMask = true(msData.numItems,1);
        end
      end
      % Prepare data objects to pass to decomposer
      if all(itemMask)
        % Full data is used
        sourceData = msData;
        if supervised
          labelData = labels;
        end
      else
        % Partial data is used
        sourceData = msData.reduce(itemMask);
        if supervised
          labelData = labels.reduce(itemMask);
        end
      end
      printFlag = false;
      if supervised
        decompositionObj = decomposer.decompose(sourceData,labelData,printFlag);
      else
        decompositionObj = decomposer.decompose(sourceData,printFlag);
      end
      % Create feature map
      varargin={'projectionType',obj.mapParams.projectionType};
      map = MSBasisMap(decompositionObj,msData.mzVector,class(obj),varargin{:});
    end
    
    function vararginParams = paramsForSpecificDecomposition(obj, params, decomposerType)
      % create cell array with optional name-value parameters
      decomposerParams=params;
      decomposerParams=rmfield(decomposerParams, 'projectionType');
      % Here extra parameters for decomposer type should be removed.
      decomposerParams=rmfield(decomposerParams, 'patterns');
      
      paramNames=fields(decomposerParams);
      nParams=length(paramNames);
      vararginParams=cell(1,2*nParams);
      for i=1:nParams
          vararginParams{2*i-1}=paramNames{i};
          vararginParams{2*i}=decomposerParams.(paramNames{i});
      end
    end
  end
end
function bool=isPosInteger(x)
bool=isPosNumber(x)&&(floor(x)==x);
end
function bool=isNonNegInteger(x)
bool=isPosInteger(x)||x==0;
end
function bool=isPosNumber(x)
bool=isscalar(x)&&isnumeric(x)&&x>0;
end
function bool=validatePattern(pattern)
    if isempty(pattern)
        bool=true;
    else
        try 
            MSPatternCell.validatePatternCell(pattern);
            bool=true;
        catch
            bool=false;
        end
    end
end
function bool=isSupervisedDecomposerType(x)
	bool=length(x)>=3&&strcmp(x(1:3),'sup');
end
