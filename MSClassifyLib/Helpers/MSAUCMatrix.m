function A = MSAUCMatrix (X, L)
  % A = AUCMatrix(X, L): Compute matrix of pairwise AUC values from score
  %   matrix X (n-by-p, columns correspond to classes) and true classes L
  %   (integer vector of length n, values 0..p, items with L==0 ignored).
  
  if ~(isnumeric(X) && ismatrix(X))
    error('Argument X must be a numeric matrix');
  end
  if ~(isnumeric(L) && isvector(L) && length(L) == size(X,1) && ...
       all(L >= 0) && all(L <= size(X,2)) && all(mod(L,1) == 0))
    error(['If argument X is n-by-p, argument L must be an integer vector ' ...
           'of length n with values in 0..p']);
  end
  
  n = size(X,2);
  A = zeros(n,n);
  for i = 1:n
    for j = [1:i-1 i+1:n]
      A(i,j) = MSROC(X(L==i,i), X(L==j,i), true)+0.5;
    end
  end
end
