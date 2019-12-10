classdef MSMassShiftProfile < matlab.mixin.Copyable
  % Mass shift analysis on MALDI peptide spectra
  %
  % Estimate mass shifts on a MALDI data object by comparing observed mass
  % defects with predictions according to the peptide mass rule.
  %
  % Properties:
  %   isValid: Validity flag
  %   mzVector: Original data m/z vector
  %   mzRange: Original data m/z range
  %   numItems: Number of profile items (e.g. for disjoint data labels)
  %   profileLength: Length of profile (number of m/z sampling points)
  %   ref: Reference m/z bin centers and corridor radius
  %   pk: Local maxima positions and shifts
  %   prf: Profile tables
  %
  % Methods:
  %   MSMassShiftProfile: Constructor
  %   computeProfile: Compute mass shift profile from MALDI data
  %   show: Generate mass defect diagram
  %   itemProfile: Return profile table for the k-th item
  %   profileColumn: Return profile table for specified column
  %   nominalToExactMass: Convert nominal to exact mass
  %   massShift: Compute mass shifts by interpolation
  %   signalRange: Compute peptide signal range
  %   setOption: Specify parameter(s) as name/value pair(s)
  %
  % MSMassShiftProfile uses the handle semantic, i.e. when assigning an
  % object of this class to a variable, only a reference to the original
  % object is copied. Use the copy method to create a (shallow) copy.
  
  properties (SetAccess = immutable)
    % Peptide mass defect rule parameters
    mdthRel = 4.95e-4;        % Relative mass defect
    mdthSig = [2.0e-5 0.02];  % Linear and constant term for stdev model
  end
  
  properties (SetAccess = protected)
    % Object status
    isValid = false; % Set to true after computing the profile
    
    % Parameters
    p = struct;  % Analysis parameters
    style = struct; % Plot style parameters
    
    % Results
    mzVector; % Original data m/z vector
    ref; % Reference m/z bin centers and corridor radius
    pk = {}; % Local maxima positions and shifts
    prf = {}; % Profile tables
  end
  
  properties (Dependent)
    mzRange; % Original data m/z range
    numItems; % Number of profile items
    profileLength; % Length of profile (number of m/z sampling points)
  end
  
  methods
    function obj = MSMassShiftProfile (varargin)
      % Constructor
      % obj = MSMassShiftProfile: Create empty object with default parameters
      % obj = MSMassShiftProfile(msData): Generate profile from msData
      %   using default parameters
      % obj = MSMassShiftProfile(__, name, value, ...): Specify parameters
      %   as name value pairs (see method setOption())
      
      % Check whether first argument is an MSMaldiData object
      msData = [];
      if nargin >= 1
        if isa(varargin{1}, 'MSMaldiData')
          % Store msData argument and remove from argument list
          msData = varargin{1};
          varargin(1) = [];
        end
      end
      
      % Initialize parameters from argument list
      obj.setOption(varargin{:});
      % Initialize plot styles
      obj.initStyles();
      
      if isfinite(obj.p.massScaleDelta)
        obj.mdthRel = obj.p.massScaleDelta;
      end
      
      % Process input data
      if ~isempty(msData)
        obj.computeProfile(msData);
      end
    end
    
    function mzRange = get.mzRange (obj)
      % Get original data m/z range
      if isempty(obj.mzVector)
        mzRange = [];
      else
        mzRange = obj.mzVector([1 end]);
      end
    end
    
    function N = get.numItems (obj)
      % Get number of profile items
      N = numel(obj.prf);
    end
    
    function N = get.profileLength (obj)
      % Get profile length
      N = size(obj.ref, 1);
    end
    
    function setOption (obj, varargin)
      % Set parameters
      % obj.setParameters(name, value, ...): Set parameters from name
      %   value pairs. The following parameters are supported:
      %
      % nNeighbors: Number of neighbors used in each direction for
      %   finding local maxima (integer, default = 0 for automatic choice)
      % massBins: Number of intervals for subdivding the mass axis
      %   (integer, default = 50)
      % fSigma: Scaling factor in terms of mutliples of standard
      %   deviation to be used for reference corridor and estimation
      %   (positive number, default = 2.0)
      % fLimits: Threshold value for ratio between estimated spread and 
      %   reference corridor. If this is exceeded, signal is considered
      %   to be noise (positive number, default = 1.2).
      % massScaleDelta: Relative mass defect shift per mass unit
      %   (default = 4.95e-4, specific to peptides)
      % unwrap: Specify whether phase unwrapping is applied (default: false)
      %
      
      % Validation functions
      isScalar = @(x) isscalar(x) && isnumeric(x);
      isScalarGT0 = @(x) isScalar(x) && x > 0;
      isScalarIntGT0 = @(x) isScalarGT0(x) && mod(x,1) == 0;
      isScalarIntGE0 = @(x) isScalar(x) && x >= 0 && mod(x,1) == 0;
      isLogicalScalar = @(x) isscalar(x) && islogical(x);
      % Setup parser
      params = inputParser;
      params.addParameter('nNeighbors', 0, isScalarIntGE0);
      params.addParameter('massBins', 50, isScalarIntGT0);
      params.addParameter('fSigma', 2, isScalarGT0);
      params.addParameter('fLimits', 1.2, isScalarGT0);
      params.addParameter('massScaleDelta', nan, isScalar);
      params.addParameter('unwrap', false, isLogicalScalar);
      
      % Parse input arguments
      params.parse(varargin{:});
      
      % Populate parameters property
      pNames = fieldnames(params.Results);
      % If options are already defined, ignore defaults
      if ~isempty(fieldnames(obj.p))
        pNames = pNames(~ismember(pNames, params.UsingDefaults));
      end
      % Store into parameters property
      for k = 1:length(pNames)
        obj.p.(pNames{k}) = params.Results.(pNames{k});
      end
    end
    
    function computeProfile (obj, msData, itemMask)
      % Compute mass shift profile
      % obj.computeProfile(msData): Compute profile (single item) for
      %   msData (MSMaldiData object)
      % obj.computeProfile(msData, itemMask): Compute separate profile
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
      itemMask = MSCheckItemMask(msData, itemMask);
      
      % Initialize result properties
      nMasks = size(itemMask,2);
      obj.mzVector = msData.mzVector(:)'; % Must be a row vector
      obj.setReference();
      obj.pk = cell(1,nMasks);
      obj.prf = cell(1,nMasks);
      
      % Iterate over mask items
      for kMask = 1:nMasks
        % Compute mean data
        if all(itemMask(:,kMask))
          % Avoid duplicating data in case all items are included
          X = mean(msData.data,1);
        else
          X = mean(msData.data(itemMask(:,kMask),:),1);
        end
        % Compute peptide shifts for local maxima
        obj.setPeptideShifts(X, kMask);
        % Compute mass shift profile
        obj.setProfile(kMask);
      end
      
      % Set valid flag
      obj.isValid = true;
    end
    
    function show (obj, K)
      % Plot mass shift profile(s)
      
      obj.checkValid;
      if nargin < 2
        K = [];
      end
      K = obj.checkItemIndex(K, true);
      
      if length(K) > 1
        % Multiple items, arrange in subplots
        tiles = MSFitTiling(length(K), 4/3);
        clf;
        for k = 1:length(K)
          subplot(tiles(1), tiles(2), k);
          obj.showProfile(k);
        end
      elseif length(K) == 1
        % Single item, plot in current axes
        obj.showProfile(1);
      end
    end
    
    function T = itemProfile (obj, k)
      % Return profile table for the k-th item (default: k = 1)

      obj.checkValid;
      if nargin < 2
        k = [];
      end
      k = obj.checkItemIndex(k, false);
      
      T = [obj.ref obj.prf{k}];
    end
    
    function X = profileColumn (obj, col, K)
      % Return profile table for specified column
      % X = obj.profileColumn(col): Return profile table with column 'col'
      %   for all items. Column may be given as column index or name
      %   (default: col = 'median')
      % X = obj.profileColumn(cok, K): Return profile table for profile
      %   items K (multiple indices allowed)
      
      obj.checkValid;
      if nargin >= 2
        if ~(ischar(col) || (isscalar(col) && isnumeric(col) && col > 0 && mod(col,1) == 0))
          error('Argument must be a string or an integer scalar')
        end
      else
        col = 'median';
      end
      if nargin < 3
        K = [];
      end
      K = obj.checkItemIndex(K, true);

      X = zeros(obj.profileLength, length(K));
      for j = 1:length(K)
        X(:,j) = obj.prf{K(j)}{:,col};
      end
    end
    
    function Y = nominalToExactMass (obj, X)
      % Convert nominal mass X to exact mass Y according to the peptide
      % mass rule
      Y = (1+obj.mdthRel)*X;
    end
    
    function Y = exactToNominalMass (obj, X)
      % Convert exact mass X to closest nominal mass Y according to the
      % peptide mass rule
      Y = round(X/(1+obj.mdthRel));
    end
    
    function dmz = massShift (obj, mz, K)
      % Compute mass shifts for given m/z values by interpolation
      % dmz = obj.massShift(mz): Compute mass shifts for m/z vector mz for
      %   all profile items
      % dmz = obj.massShift(mz, K): Compute mass shifts for m/z vector mz
      %   for profile items K (multiple indices allowed)
      
      obj.checkValid;
      if nargin < 3
        K = [];
      end
      K = obj.checkItemIndex(K, true);
      
      dmz = interp1(obj.ref.mz, obj.profileColumn('median', K), mz, 'spline');
    end  
    
    function mzBounds = signalRange (obj, condition, K)
      % Compute peptide signal range
      % mzBounds = obj.signalRange(condition, K)
      %   condition: Specify condition for signal range in case of multiple
      %     profile items. Possible values:
      %     - 'all': Signal range is where all items show signal (default)
      %     - 'any': Signal range is where at least one item shows signal
      %     - n    : Signal range is where at least n items show signal
      %   K: Vector of profile item indices (default: all items)

      % Check input arguments
      narginchk(1,3);
      obj.checkValid;
      % Item index argument
      if nargin < 3
        K = [];
      end
      K = obj.checkItemIndex(K, true);
      % Condition argument
      if nargin < 2
        condition = [];
      end
      if isempty(condition) || strcmpi(condition, 'all')
        condition = length(K);
      elseif strcmpi(condition, 'any')
        condition = 1;
      elseif ~(isnumeric(condition) && isscalar(condition) && ...
               condition > 0 && condition <= length(K) && mod(condition,1) == 0)
        error('Invalid condition argument');
      end
      
      % Identify m/z positions where peptide signal is found
      isSignal = sum(obj.profileColumn('signal', K), 2) >= condition;
      % First and last m/z positions with signal
      indBounds = [find(isSignal,1) find(isSignal,1,'last')];
      if ~isempty(indBounds)
        mzBounds = zeros(1,2);
        if indBounds(1) > 1
          mzBounds(1) = obj.ref.mz(indBounds(1)+[-1 0])'*[1;1]/2;
        else
          mzBounds(1) = obj.mzRange(1);
        end
        if indBounds(2) < obj.profileLength 
          mzBounds(2) = obj.ref.mz(indBounds(2)+[0 1])'*[1;1]/2;
        else
          mzBounds(2) = obj.mzRange(2);
        end
      else
        mzBounds = [];
      end
    end
  end
  
  methods (Access = protected)
    function checkValid (obj)
      % Error abort if object is not valid
      if ~obj.isValid
        error('No valid profile')
      end
    end
    
    function K = checkItemIndex (obj, K, allowMultiple)
      % Check item index argument
      
      if isempty(K)
        if allowMultiple
          K = 1:obj.numItems;
        else
          K = 1;
        end
      else
        if ~(isnumeric(K) && all(K > 0) && all(mod(K,1) == 0))
          error('Argument k must be positive integer')
        elseif any(K > obj.numItems)
          error('Argument k out of range')
        end
        if ~(allowMultiple || isscalar(K))
          error('Argument k must be scalar');
        end
        K = K(:)';
      end
    end
    
    function initStyles (obj)
      % Init plot style parameter structs
      colors = linspecer(3);
      obj.style.coDots = colors(1,:);
      obj.style.coReference = 'k';
      obj.style.stReference = '--';
      obj.style.coEstimation = colors(2,:);
      obj.style.stEstimation = '-';
      obj.style.coLimits = colors(3,:);
      obj.style.stLimits = '-';
    end
    
    function setReference (obj)
      % Setup reference table of m/z bin centers and corridor radius
      obj.ref = table;
      % m/z bin centers
      obj.ref.mz = ((1:obj.p.massBins)'-0.5)*(obj.mzRange*[-1 1]')/obj.p.massBins+obj.mzRange(1);
      % Bin width vector for original mass axis
      mzWidth = [obj.mzVector(2)-obj.mzVector(1), ...
                 0.5*(obj.mzVector(3:end)-obj.mzVector(1:end-2)), ...
                 obj.mzVector(end)-obj.mzVector(end-1)];
      % Corresponding std deviation
      sigmaBin = interp1(obj.mzVector, mzWidth, obj.ref.mz, 'linear')/sqrt(12);
      % Compute reference corridor
      sigmaTh = polyval(obj.mdthSig, obj.ref.mz);
      obj.ref.corridor = obj.p.fSigma*sqrt(sigmaTh.^2+sigmaBin.^2);
    end
    
    function setPeptideShifts (obj, X, k)
      % Compute peptide shifts from mean spectrum X
      
      T = table;
      % Number of neighbor positions
      n = obj.p.nNeighbors;
      if n == 0
        % Apply heuristic to select n, downsample X if necessary
        % Compute bin width vector
        d = [diff(obj.mzVector([1 2])), diff(obj.mzVector([end-1 end]))];
        if d(end)/d(1) > 5
          % Bin widths differ by factor >5, resample to new mass axis
          if d(end) < 0.1
            w = d(end)*[0.2 1];
          elseif d(1) < 0.02
            w = [0.02 0.1];
          else
            w = [d(1) min(5*d(1), 0.25)];
          end
          mzv = MSMakeMzAxis(obj.mzVector([1 end]), w);
          X = interp1(obj.mzVector, X, mzv, 'pchip');
        else
          mzv = obj.mzVector;
        end
        if length(obj.mzVector) > 1
          % Apply heuristic to select n:
          %   n*max(d) <= 1/2 AND ( n*min(d) <= 1/4 BUT at least n >=2 )
          n = floor(0.5/(max(d(2), 2*min(d(1), 0.125))));
        else
          n = 2;
        end
      end
      % Compute local maxima
      iMax = obj.localMaxima(X, n);
      T.mz = mzv(iMax)';
      % Compute mass defects
      T.shift = mod(T.mz/(1+obj.mdthRel)+0.5, 1)-0.5;
      % Store as k-th peptide shift table
      obj.pk{k} = T;
    end
    
    function setProfile (obj, k)
      % Compute shift profile from local maxima shifts

      % Compute m/z bin indices for local maxima
      mzBin = interp1(obj.ref.mz, (1:obj.p.massBins)', obj.pk{k}.mz, 'nearest', 'extrap');
      % Compute median, inter quantile range, and count
      mdMed = accumarray(mzBin, obj.pk{k}.shift, [], @median);
      mdIqr = accumarray(mzBin, obj.pk{k}.shift, [], @iqr);
      mdN = accumarray(mzBin, 1, [], @sum);
      % Compute outlier bounds assuming normal distribution and fSigma
      % multiples of std
      mdBound = obj.p.fSigma*mdIqr/(2*norminv(0.75,0,1));
      % Apply [1 1 1] smoothing to median
      mdMed = unwrap(2*pi*mdMed)/(2*pi);
      mdMed = conv(mdMed([1 1:end end]), [1 1 1]/3, 'valid');
      if ~obj.p.unwrap
        mdMed = mod(mdMed+0.5, 1)-0.5;
      end
      % Apply median plus max smoothing to outlier bounds
      mdBound = median([mdBound([1 1:end-1]), mdBound, mdBound([2:end end])], 2);
      mdBound = max([mdBound([1 1:end-1]), mdBound, mdBound([2:end end])], [], 2);
      % Signal range is where bound < fLimits * refCorridor
      mSignal = mdBound < obj.p.fLimits*obj.ref.corridor;
      % Store to profile table
      obj.prf{k} = table;
      obj.prf{k}.median = mdMed;
      obj.prf{k}.bound = mdBound;
      obj.prf{k}.signal = mSignal;
      obj.prf{k}.iqr = mdIqr;
      obj.prf{k}.N = mdN;
    end
    
    function showProfile (obj, k)
      % Plot single profile item
      
      % Plot maxima locations
      if obj.p.unwrap
        F = griddedInterpolant(obj.ref.mz, obj.prf{k}.median, 'linear', 'nearest');
        y = obj.pk{k}.shift+floor(F(obj.pk{k}.mz)+0.5);
        plot(obj.pk{k}.mz, y, '.', 'Color', obj.style.coDots);
      else
        plot(obj.pk{k}.mz, obj.pk{k}.shift, '.', 'Color', obj.style.coDots);
        set(gca, 'YLim', [-0.5 0.5]);
      end
      
      set(gca, 'defaultLineLineWidth', 1);
      hold on
      
      % Plot reference line and corridor
      plot(obj.mzRange, [0 0], ...
           'Color', obj.style.coReference, 'LineStyle', obj.style.stReference);
      plot(obj.ref.mz, obj.ref.corridor*[-1 1], ...
           'Color', obj.style.coReference, 'LineStyle', obj.style.stReference);
      
      % Plot outlier bounds and mz axis resolution
      plot(obj.ref.mz, [obj.prf{k}.median obj.prf{k}.bound]*[1 0; 1 -1; 1 1]', ...
           'Color', obj.style.coEstimation, 'LineStyle', obj.style.stEstimation);
      plot(obj.ref.mz, [obj.prf{k}.median obj.ref.corridor]*[1 -1; 1 1]', ...
           'Color', obj.style.coEstimation, 'LineStyle', obj.style.stReference);
         
      % Plot signal range limits
      signalLimit = obj.signalRange([], k);
      if ~isempty(signalLimit)
        showLimit = ((signalLimit-obj.ref.mz([1 end])').*[1 -1] > 0);
        for j = find(showLimit(:)')
          plot(signalLimit(j)*[1 1], [0.5 -0.5], ...
               'Color', obj.style.coLimits, 'LineStyle', obj.style.stLimits);
        end
      end
      hold off
    end
  end
  
  methods (Static)
    function I = localMaxima (X, m)
      % Compute indices of local maxima of row vector X, looking at m
      % neighbors left and right
      assert(isrow(X));
      M = nan(size(X));
      for k = [-m:-1 1:m]
        M = max(M, [nan(1,-k) X(max(1,1+k):min(end,end+k)) nan(1,k)]);
      end
      I = find(X > M);
    end
    
  end
  
end

