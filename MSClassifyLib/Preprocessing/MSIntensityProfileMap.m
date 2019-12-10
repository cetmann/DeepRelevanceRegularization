classdef MSIntensityProfileMap < matlab.mixin.Copyable
  % Generate intensity profile of MALDI spectra
  %
  % Estimate intensity distributions of MALDI datasets on different
  % intervals of the full m/z range
  %
  % Methods:
  %   MSIntensityProfileMap: Constructor
  %   map: Compute intensity profile from MALDI data
  %   setOption: Specify parameter(s) as name/value pair(s)
  %
  % The class uses the handle semantic, i.e. when assigning an object of
  % this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a (shallow!) copy.
  
  properties (SetAccess = protected)
    % Parameters
    p = struct;  % Analysis parameters
    
    % Profile data properties
    mzCenters; % Vector of m/z interval centers
    mzBounds;  % Vector of m/z interval bounds
    qScale;    % Quantile scale vector
    qData;     % Matrix of scaling factors to reproduce data quantiles
  end
  
  properties (Dependent)
    nMz; % Number of m/z intervals
    nQ;  % Number of points on quantile scale
  end
  
  methods
    function obj = MSIntensityProfileMap (varargin)
      % Constructor
      % obj = MSIntensityProfileMap: Create object with default parameters
      % obj = MSIntensityProfileMap(name, value, ...): Specify parameters
      %   as name value pairs (see method setOption())
      
      % Initialize parameters from argument list
      obj.setOption(varargin{:});
    end
    
    function nMz = get.nMz(obj)
      % Get number of m/z intervals
      nMz = length(obj.mzCenters);
    end
    
    function nQ = get.nQ(obj)
      % Get number of points on quantile scale
      nQ = length(obj.qScale);
    end
    
    function setOption (obj, varargin)
      % Set parameters
      % obj.setParameters(name, value, ...): Set parameters from name
      %   value pairs. The following parameters are supported:
      %
      % mzRanges: Specify m/z ranges
      %   - positive integer: number of m/z intervals (default = 5)
      %   - vector: interval bounds (positive, strictly increasing)
      % qScale: Specify quantile scale
      %   - empty: sigmoid scale computed from the data
      %   - positive integer: sigmoid scale with given number of sampling points
      %   - vector: explicit quantile scale (strictly increasing in 0..1)
      
      % Validation functions
      isScalarInt = @(x) isscalar(x) && isnumeric(x) && x > 0 && mod(x,1) == 0;
      isPosIncreasing = @(x) isvector(x) && isnumeric(x) && all(x >= 0) && ...
                             all(diff(x) > 0);
      validQScale = @(x) isempty(x) || isScalarInt(x) || isPosIncreasing(x);
      validMzRanges = @(x) isScalarInt(x) || isPosIncreasing(x);
      % Setup parser
      params = inputParser;
      params.addParameter('qScale', [], validQScale);
      params.addParameter('mzRanges', 5, validMzRanges);
      
      % Parse input arguments
      params.parse(varargin{:});
      
      % Populate parameters property
      obj.p = obj.setParams(params, obj.p);
      % If mzRanges/qScale is a vector, force to row vector
      if isvector(obj.p.mzRanges)
        obj.p.mzRanges = obj.p.mzRanges(:)';
      end
      if isvector(obj.p.qScale)
        obj.p.qScale = obj.p.qScale(:)';
      end
    end
    
    function P = map (obj, msData, itemMask)
      % Compute intensity profile
      % P = obj.map(msData): Compute profile for MSMaldiData object
      % P = obj.map(msData, itemMask): Compute separate profile
      %   items for data subsets specified by itemMask, which may be an
      %   MSLabelData object, a numerical label vector or a logical vector
      %   or matrix.
      
      % Check input arguments
      narginchk(2,3)
      if ~isa(msData, 'MSMaldiData')
        error('Data argument must be an MSMaldiData object');
      end
      if nargin < 3
        itemMask = [];
      end
      
      % Check validity of itemMask and convert to logical
      itemMask = any(MSCheckItemMask(msData, itemMask),2);
      
      % Setup quantile scale and m/z ranges
      obj.initMzRanges(msData);
      indIntervals = obj.getIndIntervals(obj.mzBounds, msData.mzVector);
      obj.initQScale(msData, indIntervals);

      % Compute single spectra quantile functions
      [A,obj.qData] = obj.computeQuantiles(msData, itemMask, indIntervals);
      
      % Create feature data object
      P = MSFeatureData(A, 'MSIntensityProfileMap', msData, itemMask);
      % Copy properties
      P.featureInfo.intensityProfile = struct;
      prp = properties(obj);
      for k = 1:length(prp)
        P.featureInfo.intensityProfile.(prp{k}) = obj.(prp{k});
      end
      
      % Clear transient properties
      obj.clear();
    end    
  end
  
  methods (Static)
    function X = sigmoidScale (da, db, n)
      % X = sigmoidScale(da,db,n): Generate a 1-by-n vector X representing
      %   a sigmoid scale on [0..1] with spacing da and db at lower and
      %   upper end of the scale, resp. If n is omitted, it is chosen as
      %   the maximum value such that increments are still increasing/
      %   decreasing at lower/upper end.

      % Definition of helper functions
      logit = @(x) log(x./(1-x));
      sigmoid = @(x) 1./(1+exp(-x));
      % If n is not specified, compute default value
      if nargin < 3 || isempty(n)
        alMin = max(logit(2*da)-logit(da), logit(2*db)-logit(db));
        n = floor(-(logit(da)+logit(db))/alMin)+1;
        if ~(isreal(n) && isfinite(n) && n >= 3)
          n = 3;
        end
      end
      % Compute scaling parameter
      al = -(logit(da)+logit(db))/(n-1);
      % Compute sigmoid vector
      t = (0:n-1)+logit(da)/al;
      X = sigmoid(al*t);
    end
    
    function indIntervals = getIndIntervals (mzBounds, mzVector)
      % Compute index intervals for m/z interval bounds
      indIntervals = ...
        [interp1(mzVector(:), 1:length(mzVector), mzBounds(1:end-1), 'next', 'extrap'); ...
         interp1(mzVector(:), 1:length(mzVector), mzBounds(2:end), 'previous', 'extrap')];
    end
  end
  
  methods (Access = protected)
    function clear (obj)
      % Clear transient properties
      obj.mzCenters = [];
      obj.mzBounds = [];
      obj.qScale = [];
      obj.qData = [];
    end
    
    function initMzRanges (obj, msData)
      % Setup m/z vector and ranges
      mzVector = msData.mzVector(:)';  % Must be a row vector
      if isscalar(obj.p.mzRanges)
        % Regular subdivision
        dmz = diff(mzVector([1 end]))/obj.p.mzRanges;
        obj.mzBounds = (0:obj.p.mzRanges)*dmz+msData.mzVector(1);
      else
        % Explicit interval bounds
        obj.mzBounds = obj.p.mzRanges;
      end
      obj.mzCenters = (obj.mzBounds(1:end-1)+obj.mzBounds(2:end))/2;
    end
    
    function initQScale (obj, msData, indIntervals)
      % Setup quantile scale
      if isempty(obj.p.qScale) || isscalar(obj.p.qScale)
        % Sigmoid quantile scale, range adapted to data
        d0 = 0;
        d1 = 0;
        minPoints = 5;  % Min number of data points within one q interval
        % Iterate over m/z intervals
        for k = 1:size(indIntervals,2)
          ind = indIntervals(:,k);
          n = ind(2)-ind(1)+1;
          ndzero = max(sum(msData.data(:,ind(1):ind(2)) <= 0, 2));
          d0 = max(d0, max(ndzero+1, minPoints)/n);
          d1 = max(d1, minPoints/n);
        end
        % Compute sigmoid scale
        obj.qScale = obj.sigmoidScale(d0, d1, obj.p.qScale);
      else
        % Explicit quantile scale
        obj.qScale = obj.p.qScale;
      end
    end
    
    function [A,S] = computeQuantiles (obj, msData, itemMask, indIntervals)
      % Compute distribution quantile functions

      % Initialize quantile function array
      numItems = msData.numItems;
      A = zeros(numItems, obj.qMapN(obj.nQ), obj.nMz);
      S = zeros(numItems, obj.nMz);
      % Compute quantiles in chunks to avoid large temporary arrays
      maskedItems = find(itemMask);
      chunks = MSData.getChunks(length(maskedItems));
      for k = 1:obj.nMz
        % Quantile functions for single spectra
        ind = indIntervals(1,k):indIntervals(2,k);
        % Iterate over chunks
        for j = 1:size(chunks,1)
          items = maskedItems(chunks(j,1):chunks(j,2));
          [A(items,:,k), S(items,k)] = ...
            obj.qMap(quantile(msData.data(items, ind), obj.qScale, 2));
        end
      end
      A = reshape(A, numItems, []);
    end
  end
  
  methods (Static)
    function [X,S] = qMap(A)
      % Feature map applied to quantile values A
      % Output S receives scaling factors computed as L1-norm of A
      X = diff(log(max(A,0)), [], 2);
      X(isinf(X(:))) = nan;
      S = sum(abs(A),2);
    end
    
    function NY = qMapN(NX)
      % Length of quantile feature map
      NY = NX-1;
    end
    
    function A = qMapInv (X,S)
      % Inverse feature map to reproduce quantile values A
      [sx1,sx2,sx3] = size(X);
      assert(all(size(S) == [sx1,sx3]));
      A = cat(2, zeros([sx1 1 sx3]), X);
      A(isinf(A)) = nan;
      A = exp(cumsum(A, 2, 'omitnan'));
      indX0 = find(~isfinite(X));
      [i,j,k] = ind2sub([sx1 sx2 sx3], indX0);
      indA0 = sub2ind(size(A), i, j, k);
      A(indA0) = 0;
      A = A.*repmat(reshape(S, [sx1 1 sx3])./sum(abs(A),2), [1 size(A,2) 1]);
    end
 end
  
  methods (Static, Access = protected)
    function P = setParams (parser, P)
      % Copy parameters from parser (inputParser object) to struct P
      % If a non-empty struct P is passed as a second, optional input 
      % argument, only non-default parameters from parser are copied to P.

      % P defaults to an empty struct
      if nargin < 2
        P = [];
      end
      if isempty(P)
        P = struct;
      end
      
      % Obtain names of parameters stored in parser
      pNames = fieldnames(parser.Results);
      % If P is specified and non-empty, ignore defaults
      if ~isempty(fieldnames(P))
        pNames = pNames(~ismember(pNames, parser.UsingDefaults));
      end
      % Copy parameters to P
      for k = 1:length(pNames)
        P.(pNames{k}) = parser.Results.(pNames{k});
      end
    end
  end
end
