function [ROC, AUC] = MSComputeROCCurve (score, positive, weight)
  % Compute ROC curve from a list of score values
  % ROC = MSComputeROCCurve(score, positive)
  %   score:    Vector or matrix of score values
  %   positive: Logical vector identifying score entries that belong to the
  %             'positive' class
  %   weight:   Vector representing relative score weights (optional)
  %
  %   If score is a vector, length N of score and positive arguments must
  %   match. The output ROC is a N+1-by-4 matrix with the (x,y) coordinates
  %   of the resulting ROC curve (columns 1:2) and corresponding lower and
  %   upper score bounds for each point (columns 3:4).
  %   If score is a matrix, its number of rows N must match the length of
  %   the positive vector. The output ROC is a N+1-by-4-by-K array, where K
  %   is the number of score columns.
  %
  % [ROC,AUC] = MSComputeROCCurve(score, positive)
  %   In addition to the ROC curve, compute the area under curve (AUC).
  
  % Check input arguments
  if ~(isnumeric(score) && (isvector(score) || ismatrix(score)))
    error('score argument must be a numeric vector or matrix');
  end
  if isvector(score)
    score = score(:);  % Force to column vector
  end
  if ~(islogical(positive) && isvector(positive))
    error('positive argument must be a logical vector')
  elseif length(positive) ~= size(score,1)
    error('Length of positive argument must match number of rows in score');
  end
  if nargin < 3 || isempty(weight)
    weight = ones(size(score,1),1);
  elseif ~(isnumeric(weight) && isvector(weight) && all(weight >= 0))
    error('weight argument must be a non-negative vector')
  elseif length(weight) ~= size(score,1)
    error('Length of weight argument must match number of rows in score');
  elseif sum(weight) <= 0
    error('Sum of weights must not be zero');
  end
      
  sz = size(score);
  N = sz(1);
  ROC = zeros(N+1, 4, sz(2));
  AUC = zeros(1, sz(2));

  weight = weight(:)*(N/sum(weight));
  positive = positive(:);
  Npos = sum(positive.*weight);
  Nneg = N-Npos;
  for k = 1:sz(2)
    % Combine score with positive class label and sort by percentage
    sortScore = sortrows([score(:,k), positive, weight], 1);
    % Cumulative count of weights and positive diagnoses
    sortCumWgt = cumsum(sortScore(:,3));
    sortCumPos = cumsum(sortScore(:,2).*sortScore(:,3));
    % Sensitivity and specificity
    sens = (sortCumPos(end)-sortCumPos)/Npos;
    spec = (sortCumWgt-sortCumPos)/Nneg;
    % (x,y) points on ROC curve, extended by (1/1)
    % (In descending order at this point)
    xyPos = flipud([1 1; [1-spec sens]]);
    % Score intervals associated with ROC curve points
    scoreLims = flipud([[-inf; sortScore(:,1)] [sortScore(:,1); inf]]);
    % Store ROC curve points and score intervals
    ROC(:,:,k) = [xyPos scoreLims];
    if nargout >= 2
      % Area under ROC curve
      AUC(1,k) = trapz(xyPos(:,1), xyPos(:,2));
    end
  end
end
