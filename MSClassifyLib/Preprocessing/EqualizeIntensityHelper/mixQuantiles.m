function B = mixQuantiles (Q, A, P)
  % Compute quantile function for a mixed distribution
  % B = mixQuantiles(Q, A, P)
  %   Q: Vector of quantile values corresponding to the rows of A
  %   A: Matrix of column vectors representing quantile functions for 
  %      source distributions
  %   P: Vector of weights corresponding to the columns of A

  % Set tolerance threshold
  eps = 1e-10;
  
  % Check arguments
  if nargin < 3
    P = ones(1, size(A,2));
  end
  checkVector = @(x,v) assert(isnumeric(x) && isvector(x), ...
                              '%s must be a numerical vector', v);
  checkVector(Q, 'Q');
  checkVector(P, 'P');
  Q = Q(:);
  P = P(:)';
  
  assert(size(A,1) == length(Q), 'Length of Q must match number of rows of A');
  assert(size(A,2) == length(P), 'Length of P must match number of columns of A');
  assert(all(min(diff(A,1,1),[],1) > -eps), 'Columns of A must be monotonically increasing');
  assert(all(diff(Q) > 0), 'Q must be strictly monotonically increasing');
  assert(all(P >= 0), 'Entries of P must be non-negative');
  assert(sum(P) > eps, 'P must not be all zero');
  
  % Compute union of quantile values as value range of mixed distribution
  X = unique(A(:));
  % Compute cdfs of source quantile functions
  F = nan(length(X), size(A,2));
  for k = 1:size(A,2)
    F(:,k) = quantileToCDF(Q, A(:,k), X);
  end
  
  % Compute mixture of cdfs
  Y = F*P'/sum(P);
  % Compute quantile function from cdf
  B = interp1(Y, X, Q, 'linear');
  
  % Extrapolate at both ends with min and max, resp.
  extrapolate = isnan(B);
  B(extrapolate & Q < Y(1)) = X(1);
  B(extrapolate & Q > Y(end)) = X(end);
end

