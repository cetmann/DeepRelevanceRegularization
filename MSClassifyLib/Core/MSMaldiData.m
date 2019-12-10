classdef MSMaldiData < MSData
  % Store mass spec datasets
  %
  % This class stores measured mass spec datasets. It is derived from the
  % base class MSData and inherits properties for storing data items
  % (spectra) and related information (mz-vector, annotations, position
  % information).
  %
  % Properties (in addition to superclass):
  %   mzVector: Row vector of mz-values
  %   mzResolution: Mean resolution of mz-values
  %   dataNorm: Data item norms for different normalizations
  %   normalization: Current normalization type
  %   normalizationTypes: List of available normalizations
  %   metaInfo: Struct for arbitrary meta information
  %
  % Methods:
  %   MSMaldiData: Constructor
  %   setMzVector: Set mz-vector and update mzResolution (overloaded)
  %   setNormalization: Set normalization type and rescale data accordingly
  %   info: Get status information
  %   plotMean: Plot mean spectra
  %   assert: Check object consistency, abort if inconsistent
  %
  % MSMaldiData uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  
  properties (SetAccess = protected)
    mzVector = [];      % Row vector of mz-values
    dataNorm = struct;  % Data item norms for different normalizations
    normalization = ''; % Current normalization type
  end
  properties (Dependent)
    mzResolution       % Mean resolution of mz-values
    normalizationTypes % List of available normalizations
  end
  properties
    metaInfo = struct; % Struct for arbitrary meta information
  end
  
  methods
    function obj = MSMaldiData (varargin)
      % Constructor
      % obj = MSMaldiData(data, mzVector): Create data object from spectra
      %   matrix (spectra stored as rows) and mz-vector.
      % obj = MSMaldiData(slFilename): Create data object from given
      %   SCiLS Lab file. Spectra, position and region information are
      %   imported from file.
      % obj = MSMaldiData(slDumpObject): Create data object from slDump
      %   object, including spectra, position and region information.
      % obj = MSMaldiData(__, name, value, ...): Specify additional options
      %   as name value pairs:
      %   - itemMask: Logical vector indicating data items to which the 
      %               generated data object is restricted.
      %   - downsample: Integer value for position grid downsampling, i.e.
      %                 loading only every n-th value (default = 1)
      %   - mzRange: m/z range, specified as [min,max], to which the
      %              generated data object is restricted.
      %   - verbose: Verbose flag. If set to true, information messages
      %              reflecting the initialization process are displayed.
      %              Default is true for initialization from SCiLS Lab
      %              file, false otherwise.
      %   - normalization: Flag specifying whether normalizations shall be
      %                    computed (default = true).
      %   - compact: Call compact(value) on the position grid after loading
      %              the data. Specify value = false to suppress compact().
      %   - dummy: Dummy flag. If set to true, only annotations and
      %            positions are loaded. Data matrix is initialized to a
      %            single column of zeros. mzRange and normalization
      %            options are ignored.
      %
      % TODO: With verbose=false, initialization from SL file still
      % displays several messages.

      % Check number of input arguments
      narginchk(1,inf);
      % Initialize empty superclass
      obj@MSData;

      % If first argument is numeric, initialization from data matrix and
      % m/z vector is assumed
      initFromData = isnumeric(varargin{1});
      if initFromData
        argsNeeded = 2;
      else
        argsNeeded = 1;
      end
      
      % Parse optional parameters
      params = inputParser;
      isLogicalScalar = @(x) islogical(x) && isscalar(x);
      isLogicalVector = @(x) (islogical(x) && isvector(x)) || isempty(x);
      isIntScalar = @(x) isnumeric(x) && isscalar(x);
      isRange = @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 2);
      isIntVector012 = @(x) isnumeric(x) && numel(x) <= 2 && all(x >= 0) && all(mod(x,1) == 0);
      isStringPack = @(x) ischar(x) && isvector(x) && strcmp(x, 'pack');
      isCompactParam = @(x) isLogicalScalar(x) || isIntVector012(x) || isStringPack(x);
      params.addParameter('itemMask', [], isLogicalVector);
      params.addParameter('downsample', 1, isIntScalar);
      params.addParameter('mzRange', [], isRange);
      params.addParameter('verbose', ~initFromData, isLogicalScalar);
      params.addParameter('normalization', true, isLogicalScalar);
      params.addParameter('compact', false, isCompactParam);
      params.addParameter('dummy', false, isLogicalScalar);
      params.parse(varargin{argsNeeded+1:end});
      
      if initFromData
        % Initialize from data matrix and m/z vector
        if ~ismatrix(varargin{1})
          error('Data argument must be a matrix');
        elseif nargin < 2
          error('mzVector argument required when initializing from data');
        elseif params.Results.downsample > 1
          error('Downsampling not supported when initializing from data matrix');
        elseif params.Results.dummy
          error('Dummy flag not supported when initializing from data matrix');
        end
        obj.initFromData(varargin{1}, varargin{2}, params.Results);
      
      elseif ischar(varargin{1})
        % Initialize from SCiLS Lab file
        obj.initFromSLFile(varargin{1}, params.Results);
      
      elseif isa(varargin{1}, 'slDump')
        % Initialize from slDump object
        obj.initFromSLDump(varargin{1}, params.Results);
        
      else
        error('Invalid argument');
      end
      
      % Initialize normalizations
      normalize = params.Results.normalization && ~params.Results.dummy;
      if normalize
        obj.message(params.Results.verbose, 'Calculating normalizations ...');
      end
      obj.initNormalization(normalize);
      
      % Compact position grid
      isLogicEq = @(x,y) islogical(x) && x == y;
      if ~isempty(obj.positions) && ~isLogicEq(params.Results.compact, false)
        if isLogicEq(params.Results.compact, true)
          obj.positions.compact();
        else
          obj.positions.compact(params.Results.compact);
        end
      end
      
      obj.message(params.Results.verbose, 'Done');
    end
    
    function initFrom (obj, maldiData)
      % Initialize object as a full copy of maldiData
      
      if ~isa(maldiData, 'MSMaldiData')
        error('Argument must be an MSMaldiData object');
      end
      
      obj.data = maldiData.data;
      obj.annotations = maldiData.annotations;
      obj.positions = maldiData.positions;
      obj.mzVector = maldiData.mzVector;
      obj.dataNorm = maldiData.dataNorm;
      obj.normalization = maldiData.normalization;
    end
    
    function R = get.mzResolution (obj)
      % Get mean mz resolution
      if length(obj.mzVector) >= 2
        R = (obj.mzVector(end)-obj.mzVector(1))/(length(obj.mzVector)-1);
      else
        R = NaN;
      end
    end
    
    function T = get.normalizationTypes (obj)
      % Get current normalization type
      T = fieldnames(obj.dataNorm);
    end
    
    function setMzVector (obj, V)
      % Set mz-vector
      % obj.setMzVector(V): Set mz-vector. Length of vector V must match
      %   the number of columns in data matrix (dataLength).
      
      if ~(isnumeric(V) && isvector(V) && length(V) == obj.dataLength)
        error(['mzVector must be a numerical vector of length equal to the ' ...
               'number of data columns']);
      else
        obj.mzVector = V(:)';
      end
    end
    
    function setNormalization (obj, T)
      % Set normalization type and rescale data accordingly
      % obj.setNormalization(T): Set normalization type to T and rescale
      %   data matrix accordingly. T must be a string matching one of the
      %   available normalization types ('raw', 'tic', 'median', 'max').
      
      % Specified normalization type valid?
      if ~isfield(obj.dataNorm, T)
        error('Invalid normalization type, use one of (%s)', ...
              strjoin(obj.normalizationTypes, ', '));
      elseif ~strcmp(T, obj.normalization)
        % New and current normalizations are different, rescale data
        scaleFactor = obj.dataNorm.(obj.normalization)./ obj.dataNorm.(T);
        % Perform rescaling column-wise (fastest and memory efficient)
        for k = 1:obj.dataLength
          obj.data(:,k) = obj.data(:,k) .* scaleFactor;
        end
        % Update current normalization type
        obj.normalization = T;
      end
    end
    
    function varargout = info (obj)
      % Get status information
      % obj.info: Print status information
      % C = obj.info: Return status information as cell string array
      
      nargoutchk(0,1);
      % Collect status strings as cell array
      info = cell(0,2);
      info(end+1,:) = {'Data matrix', ...
                       sprintf('%d items, length = %d, total values = %d', ...
                       obj.numItems, obj.dataLength, numel(obj.data))};
      info(end+1,:) = {'Normalization', obj.normalization};
      info(end+1,:) = {'mz-vector', []};
      if ~isempty(obj.mzVector)
        info{end,2} = sprintf(['Range [%.3f, %.3f], bin width [%.4f, %.4f], ' ...
                               'mean resolution = %.4f'], ...
                        min(obj.mzVector), max(obj.mzVector), ...
                        quantile(diff(obj.mzVector), [0 1]), obj.mzResolution);
      end
      info(end+1,:) = {'Annotations', []};
      if ~isempty(obj.annotations)
        info{end,2} = sprintf('%d annotations specified', ...
                              obj.annotations.numAnnotations);
      end
      info(end+1,:) = {'Position grid', []};
      if ~isempty(obj.positions)
        info{end,2} = sprintf('%d-by-%d, %.2f%% covered', ...
                              obj.positions.gridSize, ...
                              100*obj.numItems/prod(obj.positions.gridSize));
      end
      
      if nargout == 1
        % Output argument specified, return info cell array
        varargout{1} = info;
      else
        % No output argument specified, print info
        % Find width of first column
        w = max(cellfun(@(x) length(x), info(:,1)));
        % Print info
        for k = 1:size(info,1)
          fprintf('%-*s %s\n', w+1, [info{k,1},':'], info{k,2});
        end
      end
    end
    
    function plotMean (obj, varargin)
      % Plot mean spectra
      % obj.plotMean: Plot mean spectrum of all data items
      % obj.plotMean(A): Plot means of data items included in specified 
      %   annotations, given as a vector of index values or as a regular 
      %   expression name pattern
      % obj.plotMean(L): For MSLabelData object L, plot means of
      %   data items grouped by non-zero label values in L
      % obj.plotMean(M): For logical item mask M (vector or matrix), plot
      %   means of data items corresponding to columns of M
      
      % Assert consistency
      obj.assert;
      narginchk(1,2)
      X = obj.meanData(varargin{:}).';
      if size(X,1) == 0
        warning('No data rows to plot');
      else
        figure
        plot(obj.mzVector,X)
      end
    end
    
    function varargout = reduce (obj, S)
      % Reduce data items to specified item subset
      % obj.reduce(S): Reduce data items to the subset M specified by S.
      %   S may be one of the following:
      %     - a boolean vector of length == numItems; is copied to M
      %     - a regular expression name pattern or an index vector selecting
      %       a set of annotations; M is the union of all annotations
      %     - an MSLabelData object; M is the set of all items with non-zero
      %       label values
      %   If obj includes an annotation set and/or a position grid, these
      %   are reduced consistently.
      % D = obj.reduce(S): Same as above, but create and return a new data
      %   object of the same type as obj. The object obj itself is not
      %   changed.
      % Overloads the reduce method in MSData to reduce the dataNorm
      nargoutchk(0,1)
      % Convert item selector to boolean vector
      mask = obj.resolveItemSelector(S);
      
      [varargout{1:nargout}] = reduce@MSData(obj,mask);
      
      if nargout<1
           %For each norm type reduce the data
          fields = fieldnames(obj.dataNorm);
          for i=1:numel(fields)
            obj.dataNorm.(fields{i})= obj.dataNorm.(fields{i})(mask);
          end  

      end
    end
    
    function reduceMzRange(obj,newMzRange)
        % reduce mz range of data to newMzRange = [NewMzMin NewMzMax] 
        % (check for extreme cases) or to an interval part covering only
        % cProp = newMzRange of the current mzRange (in case
        % length(newMzRange) == 1) [20% will be removed from the left, 80%
        % from the right]
        %
        % the obj.data and obj.mzVector properties will be overwritten 
        % when this method is applied!
        
        % Assert consistency
        obj.assert;
        covMin = 0.001; %minimum coverage
        covMax = 1; %maximum coverage
        %check input argument
        if any([~isnumeric(newMzRange) ~isvector(newMzRange) ...
                length(newMzRange)>2])
            error(['newMzRange either has to be a vector with 2 entries ',...
                'specifying the mz range or one number between %1.4f and %1.4f ',...
                'specifying the part covered by the new mz range compared to',...
                'the old mz range'],covMin,covMax);
        end
        oldMzRange = [min(obj.mzVector) max(obj.mzVector)];
        weightRemLeft = 0.20;
        weightRemRight = 1-weightRemLeft;
        wProp = weightRemLeft/weightRemRight;
        if length(newMzRange) == 1
            cProp = newMzRange;
            if cProp < covMin
                warning('minimal proportion of new to old mz range set to %1.4f',covMin);
                cProp = covMin;
            elseif cProp >= covMax
                warning('proportion was set to >= %1.4f - the mz range remains unchanged',covMax);
                return;
            end
            newMzRangeMax = (cProp*diff(oldMzRange)+wProp*oldMzRange(2)+...
                oldMzRange(1))/(1+wProp);
            newMzRangeMin = newMzRangeMax-cProp*diff(oldMzRange);
            newMzRange = [newMzRangeMin newMzRangeMax];
        end
        indNew = obj.mzVector >= newMzRange(1) & obj.mzVector <= newMzRange(2);
        if ~any(indNew)
            error('restriction choice omitted since it would result in empty data');
        else
            newLims = [obj.mzVector(find(indNew,1)) obj.mzVector(find(indNew,1,'last'))];
            obj.data = obj.data(:,indNew);
            obj.mzVector = obj.mzVector(indNew);
            fprintf('mz range reduced to [ %5.2f    %5.2f ]  (using %3.2f %% of previous bins)\n',...
                newLims(1), newLims(2),100*sum(indNew)/numel(indNew));
        end
    end
    
    function resample(obj,data2)
      %uses MSResample (and consequently msresample from the bioinfomartics
      %toolbox) to resample the data onto the new mz vector in data2
      %
      %data2 can be a MaldiData object or a numeric vector (mz Vector)
      %specifying the (new) output mz-bins
      %
      %for simply reducing the mz-range (e.g. for generating small test
      %objects) using the method reduceMzRange instead of resample is
      %strongly recommended!

      % Assert consistency
      obj.assert;
      %check input argument
      if all([~isa(data2,'MSMaldiData'), ~(isnumeric(data2) & isvector(data2))])
          error('input argument has to be a MaldiData object or a numeric vector');
      end
      if isa(data2,'MSMaldiData')
          mzVecNew = data2.mzVector;
      else
          mzVecNew = data2;
      end
      mzVecNew = sort(unique(mzVecNew));
      disp('Resampling data to new mz/vector...');
      obj.data = MSResampleData(obj.data,obj.mzVector,mzVecNew);
      obj.mzVector = mzVecNew;
      disp('Done!');
    end
    
    function assert (obj)
      % Check object consistency, abort if inconsistent
      
      % Call superclass assert
      assert@MSData(obj);
      % Check member mzVector
      assert(length(obj.mzVector) == obj.dataLength, ...
             ['Assertion failed: Length of MSMaldiData.mzVector (%d) ', ...
              'must match number of data columns (%d)'], ...
             length(obj.mzVector), obj.dataLength);
      % Check member dataNorm
      T = fieldnames(obj.dataNorm);
      for k = length(T)
        assert(length(obj.dataNorm.(T{k})) == obj.numItems, ...
               ['Assertion failed: Length of MSMaldiData.dataNorm (''%s'': %d)', ...
                'must match number of data items (%d)'], ...
                T{k}, length(obj.dataNorm.(T{k})), obj.numItems);
      end
    end
    
    function initNormalization (obj, nrmFlag)
      % Initialize raw, tic, median and max normalizations
      % If nrmFlag is false, only raw normalization is initialized. If
      % omitted, nrmFlag is set to true.
      
      if nargin < 2
        nrmFlag = true;
      end
      
      obj.dataNorm = struct;
      % 'raw' represents unnormalized data
      obj.dataNorm.raw = ones(obj.numItems, 1, 'single');

      if nrmFlag
        % Compute normalization factors in chunks to avoid creating large 
        % temporary arrays
        chunks = obj.getChunks(obj.numItems);
        nChunks = size(chunks,1);
        
        % Setup normalization functions (expected to receive abs values)
        nrmFn = struct;
        nrmFn.tic = @(x) sum(x,2); % 'tic' = total ion count
        nrmFn.median = @(x) median(x,2);
        nrmFn.max = @(x) max(x,[],2);
        nx = obj.dataLength;
        nrmFn.rms = @(x) sqrt(1/nx*sum(x.^2,2));

        % Evalute normalization functions
        nrmTypes = fieldnames(nrmFn);
        nNrm = length(nrmTypes);
        fNrm = zeros(obj.numItems, nNrm, 'single');
        % Iterate over chunks
        for k = 1:nChunks
          items = chunks(k,1):chunks(k,2);
          % Compute absolute data values
          x = abs(obj.data(items,:));
          % Iterate over normalization types and evaluate functions
          for j = 1:nNrm
            fNrm(items,j) = nrmFn.(nrmTypes{j})(x);
          end
        end
        
        % Test for normalization factors equal zero
        fNrmZero = (fNrm == 0);
        fNrmAllZero = all(fNrmZero,1);
        if any(fNrmAllZero)
            warning(['Normalization(s) %s not available, ' ...
                     'all normalization factors are zero'], ...
                     strjoin(nrmTypes(fNrmAllZero), ', '));
        end
        
        % Fix normalization factors equal zero and store in object property
        % cases where single spectra are completely zero need to be replaced
        % by non-zero scaling factors
        for j = find(~fNrmAllZero)
          fNrm(fNrmZero(:,j),j) = 0.5*min(fNrm(~fNrmZero(:,j),j));
          obj.dataNorm.(nrmTypes{j}) = fNrm(:,j);
        end
      end
      
      % Initial normalization type is 'raw'
      obj.normalization = 'raw';
    end
    
    function showNormalization (obj, T)
      % Show normalization images
      % obj.showNormalization: Show all normalization images except 'raw'
      % obj.showNormalization(T): Show normalization image(s) specified by
      %   T, given as either a string or a cell string array.
      
      % Check input arguments
      if nargin < 2
        T = [];
      elseif ischar(T)
        T = {T};
      elseif ~iscellstr(T)
        error('Argument must be either a string or a cell string array');
      end
      if ~isempty(T)
        for k = 1:length(T)
          % Specified normalization type valid?
          if ~isfield(obj.dataNorm, T{k})
            error('Invalid normalization type (%s), use one of (%s)', ...
                  T{k}, strjoin(obj.normalizationTypes, ', '));
          end
        end
      else
        T = obj.normalizationTypes;
        T(strcmp(T, 'raw')) = [];
      end
      
      % Setup figure
      tiles = MSFitTiling(length(T), 4/3);
      figure
      colormap(obj.getColorMap);
      % Plot images
      for k = 1:length(T)
        h = subplot(tiles(1), tiles(2), k);
        obj.showImage(obj.dataNorm.(T{k}));
        h.Visible = 'off';
        title(T{k})
      end
    end
 end
  
  methods (Access = private)
    function initFromData (obj, data, mzVector, params)
      % Initialize data object from data matrix and m/z vector
      
      % Obtain item mask and m/z index range
      itemMask = obj.getItemMask(params, size(data,1));
      mzIndexRange = obj.getMzIndexRange(params, mzVector);
      % Check initialization from full data matrix
      if all(itemMask) && all(mzIndexRange == [1 size(data,2)])
        % Use unsubscripted assignment for improved performance
        obj.data = data;
        obj.setMzVector(mzVector);
      else
        % Initialization from restricted data matrix
        obj.data = data(itemMask, mzIndexRange(1):mzIndexRange(2));
        obj.setMzVector(mzVector(mzIndexRange(1):mzIndexRange(2)));
      end
    end
    
    function initFromSLFile (obj, slFilename, params)
      % Initialize data object from SCiLS Lab file
      
      % Open SCiLS Lab file
      obj.message(params.verbose, 'Open SCiLS Lab file');
      S = slDump(slFilename);
      
      % For some datasets, slDump's default heuristics yields severe
      % overlap between measurement regions. If this is the case, call
      % enlargewhitespace() to avoid overlap.
      N = S.getNumberOfSpectra();
      P = MSPositionGrid(S.getIndexMatrix());
      n = P.getMissingItems(N);
      if n > N*0.01
        warning('Detected spot overlap (%d/%d), trying to enlarge whitespace', n, N);
        S.enlargewhitespace();
      end
      
      % Initialize from slDump object
      obj.initFromSLDump(S, params);

      % Close SCiLS Lab file
      clear S;
    end
    
    function initFromSLDump (obj, S, params)
      % Initialize data object from slDump object

      % Load dataset properties
      obj.message(params.verbose, 'Loading dataset properties ...');
      % Read m/zvector, position grid and annotation set of full dataset
      mz = S.getMZvec();
      P = MSPositionGrid(S.getIndexMatrix());
      A = MSAnnotationSet(S);
      
      % Get item mask and m/z index range for data to load
      itemMask = obj.getItemMask(params, A.numItems);
      mzIndexRange = obj.getMzIndexRange(params, mz);

      % Apply downsampling to item mask
      downsampleMask = obj.getDownsampleMask(params, P.gridSize);
      if ~isempty(downsampleMask)
        itemMask(P.indexGrid(P.indexGrid > 0 & ~downsampleMask)) = false;
      end
      
      % Check for data items missing in position grid
      % Items may also be missing at the end (P.numItems < A.numItems)
      itemsCovered = [P.itemsCovered; false(A.numItems-P.numItems,1)];
      idxD = find(itemMask & ~itemsCovered);
      if ~isempty(idxD)
        % Some items are missing, issue warning and remove items from mask
        warning('%d/%d spots are missing in the index grid, data is ignored', ...
                length(idxD), A.numItems);
        itemMask = itemMask & itemsCovered;
      end
      % Keep original indices of ignored spectra in metaInfo property
      obj.metaInfo.omittedSpectraDuringImport = idxD;

      % Reduce position grid and annotation set to covered items
      if ~all(itemMask(1:P.numItems))
        [P, offset] = P.reduce(itemMask(1:P.numItems));
        if ~isempty(downsampleMask)
          downsampleMask = downsampleMask(offset(1)-1+(1:P.gridSize(1)), ...
                                          offset(2)-1+(1:P.gridSize(2)));
          P.compact(downsampleMask);
        end
      end
      if ~all(itemMask)
        A = A.reduce(itemMask);
      end
      
      if params.dummy
        % Skip loading of data, initialize data as single column of zeros
        obj.message(params.verbose, 'Loading of data skipped');
        obj.data = zeros(sum(itemMask), 1, 'single');
        obj.setMzVector(0);
      else
        % Load data, transposing from column-wise to row-wise storage
        % Loading in chunks to reduce memory needed for transpose operation
        obj.message(params.verbose, 'Loading data ...');
        nmz = mzIndexRange(2)-mzIndexRange(1)+1;
        items = find(itemMask);
        nItems = length(items);
        chunks = obj.getChunks(nItems);
        obj.data = zeros(nItems, nmz, 'single');
        for k = 1:size(chunks,1)
          kChunk = chunks(k,1):chunks(k,2);
          obj.data(kChunk,:) = S.getSpectra(items(kChunk), mzIndexRange)';
        end
        obj.setMzVector(mz(mzIndexRange(1):mzIndexRange(2)));
      end
      
      % Store properties in object
      obj.setPositions(P);
      obj.setAnnotations(A);
    end
  end
  
  methods (Access = protected)
    function showImage (obj, X)
      % Show vector X as image
      showImage@MSData(obj, X); % call superclass method
      caxis(prctile(X(:),[1 99])); % hotspot removal
    end
    
    function D = newFromTemplate (obj, data)
      % Create a new MSMaldiData object with specified data
      D = MSMaldiData(data, obj.mzVector);
      D.setNormalization(obj.normalization);
    end
  end
  
  methods (Static, Access = protected)
    function M = getItemMask (params, numItems)
      % Get itemMask parameter, defaulting to true(numItems,1)
      
      M = params.itemMask;
      if isempty(M)
        M = true(numItems,1);
      elseif length(M) ~= numItems
        error('Length of item mask must match number of spectra (=%d)', numItems);
      elseif ~any(M)
        error('Item mask is empty, no items selected');
      end
    end
    
    function M = getDownsampleMask (params, sz)
      % Get spot mask for grid downsampling (empty if no downsampling)
      
      if params.downsample > 1
        n = params.downsample;
        [x,y] = ndgrid(1:sz(1), 1:sz(2));
        M = mod(x,n) == 1 & mod(y,n) == 1;
      else
        M = [];
      end
    end
    
    function indRange = getMzIndexRange (params, mzVector)
      % Convert mzRange parameter to index range, defaulting to full range
      
      % Check mzVector argument
      isIncreasingVector = @(x) isnumeric(x) && isvector(x) && all(diff(x) >= 0);
      if ~isIncreasingVector(mzVector)
        error('m/z vector must be a numerical, monotonically increasing vector');
      end
      % Index range defaults to full range
      N = length(mzVector);
      if isempty(params.mzRange)
        indRange = [1 N];
      else
        % Largest index range included in given m/z range
        i1 = find(mzVector >= params.mzRange(1), 1, 'first');
        i2 = find(mzVector <= params.mzRange(2), 1, 'last');
        if isempty(i1) || isempty(i2)
          indRange = [];
        else
          indRange = [i1 i2];
        end
      end
    end
    
    function message (verbose, text)
      % If verbose is true, display text message
      if verbose
        disp(text);
      end
    end
  end
end
