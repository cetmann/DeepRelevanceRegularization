classdef MSDataSummary
  % Data summary class for MALDI data files
  %
  % This class is used to summarize MALDI data files by storing properties
  % of a data file, as well as derived feature data objects.
  %
  % Properties:
  %   valid: Bool flag indicating whether summary data is valid
  %   fileName: Full path of the original MALDI data file
  %   numItems: Number of data items (spots, spectra)
  %   dataLength: Length of data items (m/z bins)
  %   mzVector: m/z vector
  %   meanData: Average of all data items
  %   numAnnotations: Number of annotations
  %   normalizationData: Normalization factors as a feature data object
  %   intensityProfileData: Intensity profile data (feature data)
  %   intensityProfileComp: Intensity profile component data
  %   massShiftData: Mass shift data (feature data)
  %   massShiftComp: Mass shift component data
  %   topLevelLabel: Top level label data object
  %   topLevelMeanData: Mean data corresponding to top level label
  
  properties (SetAccess = immutable)
    fileName = '';             % Full path of the original MALDI data file
    numItems = 0;              % Number of data items (spots, spectra)
    mzVector = [];             % m/z vector
    meanData = [];             % Average of all data items
    numAnnotations = 0;        % Number of annotations
    normalizationData = [];    % Normalization factors (feature data)
    intensityProfileData = []; % Intensity profile data (feature data)
    intensityProfileComp = []; % Intensity profile component data
    massShiftData = [];        % Mass shift data (feature data)
    massShiftComp = [];        % Mass shift component data
    topLevelLabel = [];        % Top level label data object
    topLevelMeanData = [];     % Mean data corresponding to top level label
  end
  properties (Dependent)
    dataLength                 % Length of data items (m/z bins)
    mzResolution               % Average m/z bin width
    mzBinWidth                 % Minimum and maximum m/z bin width
  end
  
  methods
    function obj = MSDataSummary (fileName, varargin)
      % Constructor
      % obj = MSDataSummary(fileName): Create data summary object for MALDI
      %   data specified by fileName. If the file could not be read
      %   successfully, a warning is generated and the valid property is
      %   set to false.
      % obj = DataSummary(fileNames, name, value, ...): Specify additional
      %   options as name-value pairs:
      %   - loadOptions, meanDataOptions: Cell array with optional
      %       parameters passed to the MSMaldiData constructor and the
      %       MSData.meanData method, resp. 
      %   - intensityProfileOptions, massShiftOptions: Cell array with
      %       optional parameters passed to the MSIntensityProfileMap
      %       constructor and the MSComputeMassShiftData function, resp.
      %       Alternatively, the value may be set to 'false' to suppress
      %       computation of intensity profile and/or mass shift data.
      %   - componentDataOptions: Cell array with optional parameters
      %       passed to the MSComponentData constructor when computing
      %       intensity profile and mass shift component data.
      %   - topLevelPattern: Pattern argument passed to the
      %       getTopLevelLabel method
      %   - maxTopLevel: Maximum number of top level items accepted (0:any)
      %   - normalization: Data normalization, relevant for mean data 
      %       computation (default = 'tic')
      %   - prepFunction: Preprocessing function handle, operating on the
      %       MSMaldiData object passed as its argument. If specified, this
      %       function is called directly after loading the data.
      %   - noData: Suppress loading of spectral data, read only structural
      %       and meta data properties. If noData = true, all other options
      %       except topLevelPattern are ignored.
      
      % Check and parse input arguments
      if ~ischar(fileName)
        error('File name must be given as a character array');
      end
      params = inputParser;
      isFlag = @(x) islogical(x) && isscalar(x);
      isCellOrFlag = @(x) iscell(x) || isFlag(x);
      isString = @(x) isempty(x) || (ischar(x) && isvector(x));
      isPos0Int = @(x) isscalar(x) && isnumeric(x) && x >= 0 && mod(x,1) == 0;
      isFunc = @(x) isempty(x) || (isscalar(x) && strcmp(class(x), 'function_handle'));
      params.addParameter('loadOptions', {}, @iscell);
      params.addParameter('meanDataOptions', {}, @iscell);
      params.addParameter('intensityProfileOptions', {}, isCellOrFlag);
      params.addParameter('massShiftOptions', {}, isCellOrFlag);
      params.addParameter('componentDataOptions', {}, @iscell);
      params.addParameter('topLevelPattern', [], isString);
      params.addParameter('maxTopLevel', 0, isPos0Int);
      params.addParameter('normalization', 'tic', isString);
      params.addParameter('prepFunction', [], isFunc);
      params.addParameter('noData', false, isFlag);
      params.parse(varargin{:});

      P = params.Results;
      if islogical(P.intensityProfileOptions) && P.intensityProfileOptions == true
        P.intensityProfileOptions = {};
      end
      if islogical(P.massShiftOptions) && P.massShiftOptions == true
        P.massShiftOptions = {};
      end
      % Set default component data options, possibly superseded by
      % specified options
      P.componentDataOptions = [{'CLimMode', 'quantile', 'CLim', [0.01 0.99], ...
                                 'omitnan', 'columns', 'keepRatio', true, ...
                                 'Background', 'cmap', 'DotColor', [0 0 0]} ...
                                P.componentDataOptions(:)'];
      
      % Set file name
      obj.fileName = fileName;
      
      if ~P.noData
        % Load MALDI data
        D = MSMaldiData(fileName, P.loadOptions{:});
        % Apply prepFunction
        if ~isempty(P.prepFunction)
          P.prepFunction(D);
        end
        % Set normalization
        if ~isempty(P.normalization)
          if ismember(P.normalization, D.normalizationTypes)
            D.setNormalization(P.normalization);
          else
            warning('Normalization %s not available for MALDI data file %s', ...
                    P.normalization, fileName);
          end
        end

        % Compute summary
        obj.numItems = D.numItems;
        obj.mzVector = D.mzVector;
        obj.meanData = D.meanData(P.meanDataOptions{:});
        if ~isempty(D.annotations)
          obj.numAnnotations = D.annotations.numAnnotations;
        end
        obj.normalizationData = obj.getNormalizationData(D);
        if iscell(P.intensityProfileOptions)
          ipmap = MSIntensityProfileMap(P.intensityProfileOptions{:});
          obj.intensityProfileData = ipmap.map(D);
          obj.intensityProfileComp = ...
            MSComponentData(obj.intensityProfileData, P.componentDataOptions{:});
        end
        if iscell(P.massShiftOptions)
          obj.massShiftData = MSComputeMassShiftData(D, P.massShiftOptions{:});
          obj.massShiftComp = ...
            MSComponentData(obj.massShiftData, P.componentDataOptions{:});
        end
        obj.topLevelLabel = obj.getTopLevelLabel(D, P.topLevelPattern, P.maxTopLevel);
        if ~isempty(obj.topLevelLabel)
          obj.topLevelMeanData = D.meanData(obj.topLevelLabel, P.meanDataOptions{:});
        end
        
      else
        % noData is true, initialize via slDump object
        sl = slDump(fileName);
        A = MSAnnotationSet(sl);
        obj.numItems = A.numItems;
        obj.mzVector = sl.getMZvec()';
        obj.numAnnotations = A.numAnnotations;
        obj.topLevelLabel = obj.getTopLevelLabel(A, P.topLevelPattern, P.maxTopLevel);
        clear sl
      end
    end
    
    function N = get.dataLength (obj)
      N = length(obj.mzVector);
    end

    function X = get.mzResolution (obj)
      if length(obj.mzVector) > 1
        X = diff(obj.mzVector([1 end]))/(length(obj.mzVector)-1);
      else
        X = [];
      end
    end
    
    function X = get.mzBinWidth (obj)
      if length(obj.mzVector) > 1
        X = [diff(obj.mzVector([1 2])), diff(obj.mzVector([end-1 end]))];
      else
        X = [];
      end
    end

    function h = show (obj, varargin)
      % Show data summary graphs
      % obj.show: Show a 3-by-3 panel of data summary graphs
      % obj.show(selector): Show data summary graphs as specified by
      %   selector argument, given as a string or a cell string array.
      % obj.show(__, 'target', targetSpec): Show graph(s) in specified
      %   target: figure:  Create new figure (default)
      %           current: Plot in current figure/axes (single graph only)
      %           clear:   Clear current figure before plotting
      % h = obj.show(__): Return array of axes handles
      
      % Known graph names and default layout
      graphNames = {'meanSpectrum', 'labelImage', 'ticImage', ...
                    'massShift', 'massShiftMap', 'intensityProfileMap', ...
                    'limitMassShift', 'massShiftPlot', 'intensityProfilePlot'};

      % Check selector argument
      selector = [];
      if length(varargin) >= 1
        if iscellstr(varargin{1}) || ...
           (ischar(varargin{1}) && ismember(varargin{1}, graphNames)) || ...
           isempty(varargin{1})
          selector = varargin{1};
          varargin(1) = [];
          if ~isempty(selector) && ischar(selector)
            selector = {selector};
          end
          if ~all(ismember(selector, [graphNames, {''}]))
            error('Invalid graph selector argument');
          end
        end
      end
      if isempty(selector)
        % Show all graphs by default
        selector = graphNames;
      end
      
      % Check target argument
      params = inputParser;
      isValidTarget = @(x) ismember(x, {'figure', 'current', 'clear'});
      params.addOptional('target', 'figure', isValidTarget);
      params.parse(varargin{:});
      
      % Setup figure
      if strcmp(params.Results.target, 'current') && length(selector) > 1
        error('Cannot plot multiple graphs with target=current');
      end
      switch params.Results.target
        case 'figure'
          figure
        case 'clear'
          clf
      end
      
      if length(selector) > 1
        % Plot multiple graphs by calling show() recursively
        tiles = MSFitTiling(length(selector), 4/3);
        hAx = zeros(length(selector),1);
        for k = 1:length(selector)
          subplot(tiles(1), tiles(2), k);
          hAx(k) = obj.show(selector{k}, 'target', 'current');
        end
      else
        % Plot single graph
        showGraph = true;
        hAx = 0;
        % Set axes defaults
        legProps = {'Box', 'off', 'Location', 'best', 'Interpreter', 'none'}; 
        switch selector{1}
          case 'meanSpectrum'
            if isempty(obj.topLevelMeanData)
              plotData = obj.meanData;
            else
              plotData = obj.topLevelMeanData;
            end
            set(gca, 'ColorOrder', linspecer(size(plotData,1)), 'NextPlot', 'replaceChildren');
            plot(obj.mzVector, plotData', 'LineWidth', 1);
            xlim(obj.mzVector([1 end]));
            if ~isempty(obj.topLevelLabel)
              legend(obj.topLevelLabel.labels, legProps{:});
            end
            
          case 'labelImage'
            if isempty(obj.topLevelLabel)
              L = MSLabelData({'all'}, ones(obj.numItems,1), obj.normalizationData);
            else
              L = obj.topLevelLabel;
            end
            L.show(1, 'target', 'current');
              
          case 'ticImage'
            k = find(strcmp(obj.normalizationData.featureInfo.normalizationTypes, 'tic'),1);
            if ~isempty(k)
              obj.normalizationData.show(k, 'target', 'current');
            end
            
          case 'massShift'
            if isempty(obj.topLevelLabel)
              plotData = obj.massShiftData.meanData();
            else
              plotData = obj.massShiftData.meanData(obj.topLevelLabel);
            end
            h = MSPlotMassShiftProfile(plotData, ...
                  obj.massShiftData.featureInfo.massShiftProfile.mzVector);
            if ~isempty(obj.topLevelLabel)
              legend(h, obj.topLevelLabel.labels, legProps{:});
            end
            
          case 'massShiftMap'
            obj.massShiftComp.showMapped();
            
          case 'intensityProfileMap'
            obj.intensityProfileComp.showMapped();
            
          case 'limitMassShift'
            MSPlotMassShiftProfile(obj.massShiftComp.componentMap('limits'), ...
                                   obj.massShiftData.featureInfo.massShiftProfile.mzVector, ...
                                   obj.massShiftComp.limitColors());
            
          case 'massShiftPlot'
            obj.massShiftComp.plot(obj.topLevelLabel, 'LimitMarkerSize', true);
            
          case 'intensityProfilePlot'
            obj.intensityProfileComp.plot(obj.topLevelLabel);
            
          otherwise
            showGraph = false;
        end
        if showGraph
          hAx = gca;
        end
      end
      
      if nargout >= 1
        h = hAx;
      end
    end
  end

  methods (Static)
    function F = getNormalizationData (D)
      % Create an MSFeatureData object from the normalization factors of D
      
      nNrm = length(D.normalizationTypes);
      F = MSFeatureData(zeros(D.numItems, nNrm), 'getNormalizationData', D);
      for k = 1:nNrm
        F.data(:,k) = D.dataNorm.(D.normalizationTypes{k});
      end
      F.featureInfo.normalizationTypes = D.normalizationTypes;
    end
    
    function L = getTopLevelLabel (D, pattern, maxItems)
      % Extract top level regions/annotations
      % L = GetMainLabelData(D): Generate label data object representing
      %   the top level data partition according to the hierarchy of
      %   annotations as indicated by the '/'-separated annotation names. 
      %   If there is no such partition, L is returned empty.
      %   Otherwise, if D is an MSData object, a label data object is
      %   returned. If D is an MSAnnotationSet object, a cell string array
      %   of label names is returned.
      % L = GetMainLabelData(D, pattern): Use regular expression pattern
      %   with ()-grouping specifiers to define the annotation hierarchy.
      % L = GetMainLabelData(D, pattern, maxItems): Accept at most maxItems
      %   data items, return empty if limit is exceeded. Specify 0 to
      %   accept any number of items (default).

      % Check input arguments
      narginchk(1,3);
      if isa(D, 'MSData')
        if isempty(D.annotations)
          error('No annotations specified for data argument');
        end
        A = D.annotations;
      elseif isa(D, 'MSAnnotationSet')
        A = D;
        D = [];
      else
        error('Argument must be an MSData or MSAnnotationSet object');
      end
      
      if nargin < 2
        pattern = [];
      elseif ~(ischar(pattern) || isempty(pattern))
        error('Pattern argument must be a character array');
      end
      
      if nargin < 3
        maxItems = 0;
      elseif ~(isscalar(maxItems) && isnumeric(maxItems) && ...
               maxItems >= 0 && mod(maxItems,1) == 0)
        error('maxItems argument must be a non-negative integer scalar');
      end

      % Store all region names in cell array by hierarchy level
      regionNames = {};
      for k = 1:A.numAnnotations
        if isempty(pattern)
          % Hierarchy is defined by '/'-separators
          h = strsplit(A.annotations(k).name, '/');
        else
          % Hierarchy is defined by pattern tokens
          h = regexp(A.annotations(k).name, pattern, 'tokens', 'once');
        end
        regionNames(k,1:length(h)) = h;
      end

      % Count number of distinct regions on each hierarchy level
      maxDepth = size(regionNames,2);
      numDistinctRegions = zeros(1,maxDepth);
      for k = 1:maxDepth
        ind = ~cellfun(@isempty, regionNames(:,k));
        numDistinctRegions(k) = length(unique(regionNames(ind,k)));
      end
      % Find highest level with >1 distinct regions
      topLevel = find(numDistinctRegions > 1, 1);
      % If top level exists, generate label data object
      if ~isempty(topLevel)
        ind = find(~cellfun(@isempty, regionNames(:,topLevel)));
        [labelNames,~,sortMap] = unique(regionNames(ind, topLevel));
        labelInd = cell(length(labelNames),1);
        for k = 1:length(labelNames)
          labelInd{k} = ind(sortMap == k);
        end
        try
          if isempty(D)
            A.segmentLabels(labelInd);
            L = labelNames;
            nL = length(L);
          else
            L = D.partition(labelNames, labelInd);
            nL = L.numLabels;
          end
        catch
          L = [];
        end
      else
        % Return empty label
        L = [];
      end
      
      % Check number of label items
      if ~isempty(L) && maxItems > 0 && nL > maxItems
        L = [];
      end
    end
    
    function R = multipleDataSummary (fileNames, varargin)
      % Generate MSDataSummary objects for a list of files
      
      if ~iscellstr(fileNames)
        error('File names must be specified as a cell string array');
      end
      R = cell(size(fileNames));
      for k = 1:numel(fileNames)
        try
          R{k} = MSDataSummary(fileNames{k}, varargin{:});
        catch
          warning('MALDI data file %s could not be opened', fileNames{k});
        end
      end
    end
  end
end

