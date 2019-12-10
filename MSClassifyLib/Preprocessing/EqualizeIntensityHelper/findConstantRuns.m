function I = findConstantRuns (X, eps)
  % Find runs of constant values in data vector X
  % I = findConstantRuns(X): For numeric data vector X, return an n-by-2
  %   array I containing first and last indices of runs of constant values
  %   in X, i.e. X(I(k,1):I(k,2)) is constant.
  % I = findConstantRuns(X, eps): Consider values with difference <= eps as
  %   equal (eps must be a non-negative scalar, default = 0).
  
  if ~(isnumeric(X) && (isvector(X) || isempty(X)))
    error('Argument X must be a numeric vector');
  end
  if nargin < 2
    eps = 0;
  elseif ~(isnumeric(eps) && isscalar(eps) && eps >= 0)
    error('Argument eps must be a non-negative scalar');
  end
  
  % Find elements of X that are equal to their successor
  iD0 = find(abs(diff(X(:))) <= eps);
  % Among these, find indices that start a run
  iHead = find(diff([-inf; iD0; inf]) > 1);
  % Return start and end indices
  I = [iD0(iHead(1:end-1)) iD0(iHead(2:end)-1)+1];
end
