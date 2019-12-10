classdef MSROCCurve < matlab.mixin.Copyable
  %Class representing the ROCCurves for prediction labels. Given a certain
  %score (fraction of data classified to a given class per item) and a
  %ground truth generates and stores the corresponding ROC curves
  % Properties
  %   score:    Vector or matrix of score values.
  %   groundTruth: Logical vector identifying score entries that belong to the
  %             'groundTruth' class.
  %   ROC :   is a N+1-by-4-by-K matrix with the (x,y) coordinates
  %   of the resulting ROC curve (columns 1:2) and corresponding lower and
  %   upper score bounds for each point (columns 3:4), where K is the
  %   number of score columns.
  %   class: string indicating the class to which the ROC object corresponds
  % Methods
  %   MSROCCurves: Constructor
  properties (SetAccess=immutable)
      score; %Vector or matrix of score values.
      groundTruth;%Logical vector identifying score entries that belong to the
               %specific class.
      ROC; %  matrix with the (x,y) coordinates of the resulting ROC curve 
      class; %string indicating the class to which the ROC object corresponds
  end
  methods
      function obj=MSROCCurve(score, groundTruth, weight, classString)
          % Compute ROC curve from a list of score values
          % obj = MSComputeROCCurve(score, groundTruth, classString)
          %   score:    Vector or matrix of score values
          %   groundTruth: Logical vector identifying score entries that belong to the
          %             specific class
          %   weight: Vector representing relative score weights (optional)
          %   classString: (optional) indicates to which class the ROC
          %   analysis corresponds
          %   If score is a vector, length N of score and groundTruth arguments must
          %   match. The property ROC is a N+1-by-4 matrix with the (x,y) coordinates
          %   of the resulting ROC curve (columns 1:2) and corresponding lower and
          %   upper score bounds for each point (columns 3:4).
          %   If score is a matrix, its number of rows N must match the length of
          %   the groundTruth vector. The output ROC is a N+1-by-4-by-K array, where K
          %   is the number of score columns.
          
          % check input arguments
          narginchk(2,4)
          % For backwards compatibility, omitting the weight argument must
          % be supported. Thus, allow weight and classString arguments to
          % be in any order.
          isString = @(x) ischar(x) && size(x,1) <= 1;
          if nargin < 3
            % Initialize optional arguments empty
            weight = [];
            classString = '';
          elseif nargin == 3
            % If only one optional argument is given and is a string,
            % assume it is classString
            if isString(weight)
              classString = weight;
              weight = [];
            else
              classString = '';
            end
          elseif isString(weight)
            % If both optional arguments are given and the first is a
            % string, assume it is classString
            h = classString;
            classString = weight;
            weight = h;
          end
          if nargin<3
              obj.class='';
          else
              if ~ischar(classString)||size(classString,1)>1
                  error('The input argument <classString> must be a valid string')
              end
              obj.class=classString;
          end
          % Other input arguments checked in call of MSComputeROCCurve
          try
               obj.ROC=MSComputeROCCurve(score, groundTruth, weight);   
          catch ME
              emsg=ME.message;
              error(['There was some error creating the ROC object. ' emsg])              
          end
          obj.score=score;
          obj.groundTruth=groundTruth;
      end
      function R=reduceROC(obj,K)
        % Return reduced ROC matrix for score data K
        % R = obj.reduceROC(K): Return ROC matrix for score data K with
        %   those rows removed that correspond to 'shadow' points, i.e. 
        %   points that are clearly inferior to neighboring points.
        
        if ~(isscalar(K) && isnumeric(K) && K >= 1 && ...
             K <= size(obj.score,2) && mod(K,1) == 0)
          error('Invalid score data index K');
        end
        R = MSReduceROCCurve(obj.ROC(:,:,K));
      end
      function AUC=computeAUC(obj, reduced)
          %Computes area under the curve (AUC)
          %AUC=obj.computeAUC(reduced)
          %OUTPUT
          % AUC: list of AUC per class
          % If reduced is true or 'reduced', AUC is computed based on the
          % reduced ROC curve. Default is false.
          if nargin < 2
            reduced = false;
          elseif ischar(reduced) && strcmpi(reduced, 'reduced')
            reduced = true;
          end
          if ~(isscalar(reduced) && islogical(reduced))
            error('Invalid argument');
          end
          nScores = size(obj.score,2);
          AUC=zeros(1,nScores);
          
          for k=1:nScores
              % get (x,y) coordinates
              if reduced
                xyPos = obj.reduceROC(k);
              else
                xyPos=obj.ROC(:,:,k);
              end
              % Area under ROC curve
              AUC(1,k) = trapz(xyPos(:,1), xyPos(:,2));
          end
      end
      function optBalAcc = computeOptBalAcc(obj)
        optBalAcc = 1/2*squeeze( max( 1-obj.ROC(:,1,:)+obj.ROC(:,2,:), [], 1 ) )';
      end
      function optThr = computeOptThr(obj,fun)
        % Minimizes fun to retrieve optimal thr, default is (negative) balanced
        % accuracy with equal priors, i.e. fun = @(x)-c(1-x+y)
        if nargin < 2
          fun = @(M)( -1/2*sqrt(2)*( M(:,2)-M(:,1)+1 ) );
        end
        nScores = size( obj.score, 2 );
        optThr = zeros( 1, nScores );
        for k = 1:nScores
          Mfull = squeeze( obj.ROC(:,1:4,k) );
          M = Mfull(:,1:2);
          [~,optThr_k] = min( fun(M) );
          % Consider neighbors and score intervals to retrieve a finer thr
          % estimate within the two thr limits in Mfull(...,3:4)
          if optThr_k == 1
            optThrIntervalPos = 1;
          elseif optThr_k == size(M,1)
            optThrIntervalPos = 0;
          else
            MPart = Mfull( (-1:1)+optThr_k, 1:4);
            spans = MPart(:,4)-MPart(:,3);
            balAccVals = 1/2*(1-MPart(:,1)+MPart(:,2));
            if balAccVals(1) > balAccVals(3)
              neighborIdx = 1;
            else
              neighborIdx = 3;
            end
            diffBalAcc = balAccVals(2) - balAccVals( neighborIdx );
            multFactorAcc = max(1-20*diffBalAcc,0);
            factorSpan = spans(neighborIdx)/(100*eps+spans(2));
            multFactorSpan = 2/pi*atan(pi*factorSpan);
            if MPart(neighborIdx,4) > MPart(2,4)
              optThrIntervalPos = 0.5+0.5*multFactorAcc*multFactorSpan;
            else
              optThrIntervalPos = 0.5-0.5*multFactorAcc*multFactorSpan;
            end
          end
          intervalPosFun = @(pos,Interval)( Interval(1) + pos*(Interval(2)-Interval(1)) );
          optThr(k) = intervalPosFun( optThrIntervalPos, Mfull(optThr_k,3:4) );
        end
      end
      function S = optimalScore(obj)
          % Compute score with maximum average sensitivity
          nScores = size(obj.score,2);
          S = zeros(1,nScores);
          
          for k = 1:nScores
              % Find maximum of average sensitivity
              [~,ind] = max(mean([1-obj.ROC(:,1,k), obj.ROC(:,2,k)], 2));
              % Optimal score is center of corresponding score interval
              S(k) = mean(obj.ROC(ind,3:4,k));
          end
      end
      function plotCurves(obj,numFeatures)
          % Plots ROC curves
          % obj.plotCurves(numFeatures)
          % INPUT
          %  numFeatures
          
          %input validation
          narginchk(1,2)
          nScores = size(obj.score,2);
          if nargin < 2 || isempty(numFeatures)
              numFeatures=1:nScores;
          else
              if ~isnumeric(numFeatures)||~isvector(numFeatures)
                  error('The input <numFeatures> must be a numerical vector')
              end
              numFeatures = numFeatures(:)';
          end
          for k = numFeatures
              % Get xy positions of ROC curve, remove shadow points
              xyPos = obj.reduceROC(k);
              % Plot
              plot(xyPos(:,1), xyPos(:,2));
              hold on
          end
          hold off
          legend(cellstr(num2str(numFeatures')), 'Location', 'best', 'Box', 'off');
          title(['ROC-nfeatures (' obj.class ')']);
          xlabel('False positive rate (1-Specificity)');
          ylabel('True positive rate (Sensitivity)');
      end
  end
end