function [R, P] = MSROC (X, Y, signed)
  % Compute ROC criterion discriminating data distributions
  % R = MSROC(X,Y): For input matrices X (m-by-p) and Y (n-by-p) with same
  %   number of columns, compute the area under the ROC-curve associated 
  %   with using a threshold for discriminating between the two input data
  %   distributions. Output values R are within [0, 0.5] and represent the
  %   difference between this area and the diagonal line in the ROC graph.
  %   Values in the row vector R correspond to columns of X and Y and can
  %   be interpreted as the probabilities (minus 0.5) that values in one
  %   matrix are larger than in the other.
  % R = MSROC(X,Y,signed): If signed is true, compute the signed difference
  %   of the above area under curve value, which is within [-0.5, 0.5]. The
  %   default is signed=false.
  % [R,P] = MSROC(__): In addition to R, compute the p-values for the
  %   corresponding area under curve values. The p-value is always computed
  %   for the two-sided hypothesis, regardless of the signed argument.
  
  if nargin < 3 || isempty(signed)
    signed = false;
  end
  if ~(isnumeric(X) && isnumeric(Y) && ismatrix(X) && ismatrix(Y) && ...
       ~isempty(X) && ~isempty(Y) && size(X,2) == size(Y,2))
    error(['Both arguments must be non-empty 2D-arrays with the same ', ...
           'number of columns']);
  else
    c = size(X,2);
    m = size(X,1);
    n = size(Y,1);
    R = zeros(1,c);
    P = zeros(1,c);
    if c > 0
      % Iterate over columns
      for k = 1:c
        % Compute tie corrected ranks for X and Y
        tcRank = RankTieCorrected([X(:,k);Y(:,k)]');
        % Compute rank sum of values of X (i.e. the first m entries)
        R(k) = sum(tcRank(1:m));
      end
      % Compute ROC statistic as (q-0.5) where q is the probability that
      % an element of X is larger than an element of Y
      R = (R-m*(m+1)/2)/(m*n)-0.5;
      % If signed is false, take the absolute difference
      if ~signed
        R = abs(R);
      end
      % Compute the two sided p-value for the computed ROC statistic, using
      % the approximation that (signed) R is normal
      sigma = sqrt((m+n+1)/(12*m*n));
      P = 2*cdf('Normal',-abs(R),0,sigma);
    end
  end
end

function R = RankTieCorrected (X)
  % Compute the rank of all values in X, correcting ties by replacing
  % the ranks of multiple identical values by their median rank.
  
  % Sort values of X in ascending order and obtain sorting order
  [xSorted, xOrder] = sort(X);
  % To avoid numerical instabilities, replace eps-ties by true ties
  xEps = (xSorted(end)-xSorted(1))*eps;
  for k = 2:length(xSorted)
    if xSorted(k)-xSorted(k-1) < xEps
      xSorted(k) = xSorted(k-1);
    end
  end
  % Compute uncorrected rank by inverting the sorting order
  xRank = zeros(size(xOrder));
  xRank(xOrder) = 1:length(xOrder);
  % Compute mapping between original rank and rank within unique values
  [~, origRankMap, uniqueRankMap] = unique(xSorted);
  % Count multiple occurrences of sorted input values
  xCount = diff([origRankMap' length(X)+1]);
  % Apply tie correction to rank map by adding an offset depending on the
  % number of occurences
  correctedRankMap = origRankMap'+(xCount-1)/2;
  R = correctedRankMap(uniqueRankMap(xRank));
end
