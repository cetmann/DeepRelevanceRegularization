classdef MSBasisMap < MSFeatureMap
  % Feature map allowing different projections of basis vectors in
  % MSDecomposition objects
  %
  % Properties:
  %   decompositionObj: Decomposition object containing the basis vectors
  %
  % Methods:
  %   MSLinearBasisMap: Constructor
  %   map_impl: Compute output feature data from input data
  %
  % MSBasisMap uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a deep copy.
  %
%   properties
%       mapParam=struct;
%   end
  properties (SetAccess = protected)
    decompositionObj  % Decomposition object (contains basis)
    projectionType = 'correlation'; % If possible... compare constructor
  end
  properties(Constant)
      listProjectionTypes = {'scalarProduct', 'correlation', 'nonLinearGeneric'};
  end
  properties(Dependent)
      basis
      l1l2
      gini
      auc
      pvalues
  end  
  methods % Property get
      function value = get.basis(obj)
          value = obj.decompositionObj.basis;
      end
      function value = get.l1l2(obj)
          value = obj.decompositionObj.addInfo.l1l2Norms(obj.decompositionObj.sortIndexCurr);
      end
      function value = get.gini(obj)
          value = obj.decompositionObj.addInfo.giniIndices(obj.decompositionObj.sortIndexCurr);
      end
      function value = get.auc(obj)
          if ~isfield(obj.decompositionObj.addInfo,'rocMaxVals')
              value = [];
              disp('AUC values have not been computed')
              return
          end
          value = obj.decompositionObj.addInfo.rocMaxVals(obj.decompositionObj.sortIndexCurr);
      end
      function value = get.pvalues(obj)
          if ~isfield(obj.decompositionObj.addInfo,'pValsMax')
              value = [];
              disp('AUC values have not been computed')
              return
          end
          value = obj.decompositionObj.addInfo.pValsMax(obj.decompositionObj.sortIndexCurr);
      end
  end
  methods
    function obj = MSBasisMap (msDecomposition, mzVector, creator,varargin)
      % Constructor
      % obj = MSBasisMap(msDecomposition, mzVector, creator)
      %   msDecomposition: msDecomposition object or numeric matrix
      %     (automatically generating a dummy decomposition object)
      %   mzVector: optional, usually already defined in the
      %     msDecomposition object, if not specified or empty, i.e. [], it
      %     will use that one instead or 1:dataLength if not available. If
      %     mzVector is specified, dimensions are checked and will be saved
      %     instead
      %   creator: Feature map creator string (optional, defaults to 
      %     MSBasisMap)
      
      narginchk(1,inf);
      if nargin < 3|| isempty(creator)
        % If creator string is not specified, use this class name
        creator = 'MSBasisMap';
      end
      if nargin < 2 || isempty(mzVector)
        if ~isempty(msDecomposition.mzVector)
          mzVector = msDecomposition.mzVector;
        else
          error('An mzVector must be specified')
        end         
      else
        if ~(isnumeric(mzVector) && isvector(mzVector) && ...
            length(mzVector) == msDecomposition.dataLength)
          error('mzVector has to be a numeric vector with matching data length');
        end
      end
     
      % Initialize feature map
      obj@MSFeatureMap(mzVector, creator);
      % Check input arguments
%       if isempty(basis) || ~isnumeric(basis) || ~ismatrix(basis)
%         error('basis must be a non-empty numeric matrix');
%       end  
%       if size(basis,2)~=length(obj.mzVector)
%           error(['the number of columns of basis and the number'...
%                  'of elements of mzVector must be the same'])
%       end  

      % optional parameter-Define parameters for iteratively determined mappings
      params=inputParser;
%       params.addParameter('checkT', 5, @isPosInteger);
%       params.addParameter('minIt', 20, @isPosInteger);
%       params.addParameter('maxIt', 100,@isPosInteger);
%       params.addParameter('stopLim', 0.5*5.0*10^(-3), @isPosNumber);
      params.addParameter('projectionType','correlation',...
          @(x) ischar(x)&& any(strcmp(MSBasisMap.listProjectionTypes,x)));
      %parse actual parameters
      params.parse(varargin{:});
      %update properties
      obj.projectionType=params.Results.projectionType;      
%       obj.mapParam=params.Results;
%       obj.mapParam = rmfield(obj.mapParam,'projectionType');
      obj.decompositionObj = msDecomposition;
      %Previous checking of compatibility
      if strcmp(obj.projectionType,'nonLinearGeneric')&& checkNonLinearGenericButNotChild('nonLinearGeneric',...
          obj.decompositionObj.decomposerType)
        obj.projectionType = 'correlation'; % Otherwise nonLinearGeneric!
      end
      
      % stopping when relative Frobenius norm changes in K go below stopLim
      % ... compare MSDecomposer
    end
    
    function show(obj, type, varargin)
        cla
%         axis off
        narginchk(1,inf)
        %Extra inputs
        params=inputParser;
        params.addParameter('nFeatures', inf, @isPosInteger);
        params.addParameter('threshold', 0, @isNonNegNumber);
        params.parse(varargin{:});  
        params=params.Results;
        
        params.nFeatures=min(params.nFeatures, obj.numFeatures);
        if nargin < 2
            type = 'auc';
        elseif ~ischar(type)
            error('The input <type> must be a string indicating type of visualization')
        end
        do = obj.decompositionObj;
        info = do.addInfo; 
        index = do.sortInd.(do.sortTypeCurr);
        index = index(1:params.nFeatures);
        % Number of basis vectors
        switch type
            case {'auc', 'pvalue'}
                if ~isfield(info,'rocMaxVals')
                    error('No roc values have been computed to be shown')
                else       
                    discrimClass = info.discrimClass(index);
                    classes = unique(discrimClass);
                    if strcmp(type, 'auc')
                        fgValues = obj.auc;
                    elseif strcmp(type, 'pvalue')
                        fgValues = -log(obj.pvalues);
                    end
                    fgValues = fgValues(1:params.nFeatures);
                    % apply threshold
                    maskThreshold = abs(fgValues)>=params.threshold;
                    index = index(maskThreshold);
                    discrimClass = discrimClass(maskThreshold);
                    fgValues = fgValues(maskThreshold);
                    
                    nclasses = length(classes);
                    v = zeros(1,nclasses);
                    % plot prominent roc values (highest absolute value)
                    for i=1:nclasses
                        mask = discrimClass == classes(i);
                        v(i) = stem(index(mask),...
                            fgValues(mask),'filled','linewidth', 1); hold on
                    end
                    axis on
                    hold off
                    if nclasses > 1 || any(classes > 1)
                        legend(v, num2str(classes));
                    end
                end
            case 'basis'
                B = obj.basis(1:params.nFeatures,:);
                Bmax = max(B, [], 2);
                BNorm = -B./repmat(Bmax, [1 size(B,2)])+repmat((1:size(B,1))', [1 size(B,2)]);
                plot(do.mzVector, BNorm', 'k-');
                axis on
                axis tight ij; 
            case 'heat'
                % Number of m/z bins
                numMzBins = 200;
                % Map m/z index to bin
                mzRange = do.mzVector([1 end]);
                mzBinMap = max(ceil((do.mzVector-mzRange(1))/(mzRange(2)-mzRange(1))*numMzBins), 1);
                % Compute heat map
                nmfHeatMap = zeros(params.nFeatures, numMzBins);
                fnEnergy = @(x) norm(x,1);
                % fnNorm = @(x) norm(x,1);
                fnNorm = @max;
                for k = 1:params.nFeatures
                  x = accumarray(mzBinMap', obj.basis(k,:), [], fnEnergy)';
                  nmfHeatMap(k,:) = x/fnNorm(x);
                end
                position=get(gca,'Position');
                axis off
                axes('Position', [position(1) position(2) 15/16*position(3) position(4)]);
                imagesc(nmfHeatMap)
                colormap(viridis)                
                axes('Position', [position(1) + 15/16*position(3) position(2) 1/16*position(3) position(4)]);
                bScores = zeros(3,params.nFeatures);
                if ~isempty(obj.auc)
                    bScores(1,:) = 2*obj.auc(1:params.nFeatures);
                end
                bScores(2,:) = 1-obj.gini(1:params.nFeatures);
                bScores(3,:) = obj.l1l2(1:params.nFeatures);
                bScores(3,:) = bScores(3,:)/max(bScores(3,:));
                plot(bScores', repmat((1:params.nFeatures)', [1 3]), 'LineWidth', 1); 
                axis tight
                ax = gca;
                ax.YDir = 'reverse';
                set(ax,'yticklabel',[]);
            case 'signed'
               % divide nmf in positive and negative auc
                R = obj.auc(1:params.nFeatures);
                thres = params.threshold;
                underThres = abs(R) > thres;
                positiveROC = R > 0;
                negativeROC = R < 0;
                class1nmf = obj.basis(underThres & positiveROC, :);
                class2nmf = obj.basis(underThres & negativeROC, :);
                hold on
                rocValuesPos = R(R>0);
                legendCell = cell(size(class1nmf,1), 1);
                for i= 1:size(class1nmf,1)
                    plot(obj.mzVector, class1nmf(i, :))
                    legendCell{i} = (num2str(rocValuesPos(i)));
                end
                axis on
                %plot mean
                plot(obj.mzVector, mean(class1nmf, 1), 'k','LineWidth',2)
                %legend(legendCell)

                % Negatives
                rocValuesNeg = R(R<0);
                legendCell = cell(size(class2nmf,1), 1);
                for i= 1:size(class2nmf,1)
                    plot(obj.mzVector, -class2nmf(i, :))
                    legendCell{i} = (num2str(rocValuesNeg(i)));
                end
                %plot mean
                plot(obj.mzVector, -mean(class2nmf, 1), 'k','LineWidth',2) 
            case 'image'
                position=get(gca,'Position');
                for k = 1:4                    
                    axes('Position', [position(1)+0.2*(k-1) position(2)+0.6 1/5*position(3) 1/5*position(4)]);
                    set(gca,'XTickLabel','')
                    imagesc(do.coeffs.positions.encube(do.coeffs.data(:,k)),...
                    'alphadata',logical(do.coeffs.positions.indexGrid),...
                    'clipping','off'); colormap(jet);
                  axis image; axis off;
                  title(sprintf('Feature %2g',k));
                end
                axes('Position', [position(1) position(2) position(3) 1/3*position(4)]);
                set(gca,'XTickLabel','')
                M = do.basis(1:4,:)';

                h = plot(do.mzVector,M,'LineWidth',2);
                legend({'Feature 1','Feature 2','Feature 3','Feature 4'});
                
                % get same colors to plot the convolution
                c = cell2mat(get(h,'Color'));
                set(gca, 'ColorOrder',c,'NextPlot','replacechildren');

                hold on;
                mzResolution = (obj.mzVector(end)-obj.mzVector(1))/(length(obj.mzVector)-1);
                numConv = ceil(20/mzResolution); % 20 Da convolution
                MC = convn(M,1/numConv*ones(numConv,1),'same');

                plot(do.mzVector,-5*MC,'lineWidth',2); % Times 5 to see more clearly
                ylim([-5*prctile(MC(:),99.9) prctile(M(:),99.9)]);

                h = zoom;
                set(h,'Motion','horizontal','Enable','on');
            otherwise
                error('Unknown type of visualization')
        end
    end
    function MapVisualizationTool(obj)
        MapVisualizationControl('feature_map', obj);
    end
    
    function switchProjectionType(obj,newProjectionType)
      if ~ischar(newProjectionType)
        error('newProjectionType has to be a string');
      else
        if ~any(strcmp(obj.listProjectionTypes,newProjectionType))
          error('Unknown projection type');
        else
          if checkNonLinearGenericButNotChild(newProjectionType,...
              obj.decompositionObj.decomposerType)
            error(['The decomposerType of the decompositionObj property'...
              ' is not a class name for a child of MSDecomposer']);
          end
        end
      end
      obj.projectionType = newProjectionType;
    end
  end
  methods(Static)
      function string = visualizationDescription(type)
          % Returns an explanation of the give type of map visualization
          % INPUT
          %  type: One of the visualization types
          if  ~ischar(type)
              warning('The visualization type must be a string')
          end
          switch type
              case 'auc'
                  string = ['Signed AUC values corresponding to the mz-values. '...
                          'The color of each pin corresponds to the class that '...
                          'is discriminated from the rest using the associated mz-value.'];
              case 'pvalue'
                  string = ['p-values corresponding to the mz-values. '...
                          'The color of each pin corresponds to the class that '...
                          'is discriminated from the rest using the associated mz-value.'];
              case 'image'
                  string = ['Representation of the first 4 vector basis (lower figure)'...
                            'and their correspondent channels (upper figures). The '...
                            'lines with negative y-axes correspond to a convolution of'...
                            ' the vectors in wider mz-windows showing the general behaviour'...
                            ' of the intensity along the mz-vector'];
              case 'heat'
                  string = sprintf(['Representation of the basis vectors through heat maps. '...
                            'Each row of the matrix corresponds to a basis vector. The'...
                            ' color represents the average intensity of the mz-values '...
                            'corresponding to an mz-window. \n\nAt the right side of the heat '...
                            'matrix values related to AUC, gini and l1l2 for each basis '...                            
                            'vector are shown:\n \nblue->2*AUC\n red-> 1-gini\n orange-> l1l2']);
              case 'signed' 
                  string = ['Representation of the basis vectors. Negative(positive)'...
                            'vectors correspond to negative(positive) AUC  values. Mean values'...
                            ' of both, positive and negative vectors, are also included (bold black)'];
              case 'basis'
                  string = 'Compressed representation of the basis vectors.';
              otherwise
                  string = 'Unknown type of visualization';
                  warning('The requested type of visualization is not valid')
          end
      end
  end

  methods(Access = protected)
      
    function cp = copyElement(obj)
        % overloads the copyElement method in matlab.mixin.Copyable
        % makes a deep copy of the object (copies the decomposition). This
        % changes the behaviour of the method copy, which calls copyElement
        
        decompositionObjCopy=obj.decompositionObj.copy; % Decomposition object (contains basis)
        cp=MSBasisMap(decompositionObjCopy,obj.mzVector,obj.creator);
        cp.projectionType = obj.projectionType;
%         cp.mapParam=obj.mapParam;
    end
    
    function numFeatures = getNumFeatures (obj)
      numFeatures = size(obj.decompositionObj.basis,1);
    end
    
    function featureData = map_impl (obj, msData, itemMask, numFeatures, printFlag)
      % Compute output feature data from input data
      % featureData = obj.map(msData, itemMask, numFeatures, printFlag)
      %   msData: MSMaldiData object with input data (spectra, feature vectors)
      %   featureData: Resulting MSFeatureData object with feature vectors
      %                consisting of scalar products between input data
      %                items and the map's basis vectors
      %   itemMask: Logical mask specifying the spectra to which the mapping
      %             is applied
      %   numFeatures: Restrict the output number of features to the
      %               specified quantity
      %   printFlag: If true, print progress information
      
      % Create output feature data
      
      Y = msData.data(itemMask,:);
      Y = MSDecomposer.correctPrecision(Y,obj.decompositionObj.param.precision);
      X = obj.decompositionObj.basis(1:numFeatures,:);
      X = MSDecomposer.correctPrecision(X,obj.decompositionObj.param.precision);
      switch obj.projectionType
        case 'scalarProduct'
          K = Y*X.';
        case 'correlation'
%           l2normedRows = @(X)(X./repmat(sqrt(sum(X.*X,2)),[1 size(X,2)]));
%           K = l2normedRows(Y)*l2normedRows(X).';
            K = corr(Y.',X.');
        case 'nonLinearGeneric'
%           fh = str2func([obj.decompositionObj.decomposerType ...
%             '.updateRule_impl_K']);
          
          if printFlag
            fprintf('Initiating... ')
          end
          
          if ismethod(obj.decompositionObj.decomposer,'nonLinearGenericBasisMapInit')
            K = obj.decompositionObj.decomposer.nonLinearGenericBasisMapInit(Y,X);
          else
            %K = 10*eps+abs(Y/X); 
            K = 2*eps+abs(randn(size(Y,1),size(X,1)));
          end

          if printFlag
            fprintf('Done!\n');
          end

          K2 = K+1;
          t = 0;
          if printFlag
            fprintf('Updating...\n');
            timeStarted = clock;
          end
          
            YXt = Y*X';
            XXt = X*X';
          while ~MSDecomposer.stoppingCriterion(t,K,K2,0.1,0.1,obj.decompositionObj.param.checkT,...
              obj.decompositionObj.param.stopLim/2,obj.decompositionObj.param.minIt,obj.decompositionObj.param.maxIt) 
%           The param.maxIt should be replaced by a new parameter
%           maxItProjection in the future
            K2 = K;

            K = MSDecomposerNMFFrob.updateRule_impl_K_SpecifiedNLG(...
              obj.decompositionObj.param,YXt,XXt,Y,K,X,[]);
            K = MSDecomposer.projectOnPositiveDoubleOrSingle(K);
           
            t = t+1;
            if printFlag
              elapsed = etime(clock,timeStarted);
              if t==0+1
                nnns = 0;
              end
              fprintf([repmat('\b',1,nnns) '      %5d out of max. %5d it.,',...
                      '\n %8.2f s elapsed, up to %8.2f s remaining\n'],t,...
                      obj.decompositionObj.param.maxIt,elapsed,((obj.decompositionObj.param.maxIt)/t-1)*elapsed);
              nnns = 83; % String length ('overwriting' last printed line)
            end
          end
          
          if printFlag
            fprintf('Done!\n');
          end
      end
      % Create empty output feature data
      featureData = MSFeatureData(zeros(msData.numItems, numFeatures, 'like', msData.data), ...
                                  class(obj), msData, itemMask);
      % Store feature data for items selected by itemMask
      featureData.data(itemMask,:) = K;
      
    end
    
  end

end

% Local functions
function bool=isPosInteger(x)
bool=isPosNumber(x)&&(floor(x)==x);
end
function bool=isPosNumber(x)
bool=isscalar(x)&&isnumeric(x)&&x>0;
end
function bool=isNonNegNumber(x)
bool=isscalar(x)&&isnumeric(x)&&x>=0;
end
function forbidden = checkNonLinearGenericButNotChild(newProjectionType,child)
% Usually the child obj.decompositionObj.decomposerType will be 'checked'

forbidden = strcmp(newProjectionType,'nonLinearGeneric') && ...
              ~any(strcmp('MSDecomposer',...
              superclasses(child)));
            
end
