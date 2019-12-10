function F = quantileToCDF (Q, A, X)
  % Compute cumulative distribution function for given quantile function
  % F = quantileToCDF (Q, A, X)
  %   Q: Vector of quantile values (in [0..1])
  %   A: Corresponding vector of sample data values
  %   X: Vector of data values for which CDF is requested
  %   F: Resulting vector of CDF values
  % Q and A must have same length. A must be monotonically increasing, Q
  % and X must be strictly monotonically increasing.
  
  % Set tolerance threshold
  eps = 1e-10;
  
  % Check arguments
  checkVector = @(x,v) assert(isnumeric(x) && isvector(x), ...
                              '%s must be a numerical vector', v);
  checkVector(Q, 'Q');
  checkVector(A, 'A');
  checkVector(X, 'X');
  assert(length(Q) == length(A), 'Lengths of Q and A must match');
  assert(all(Q >= 0 & Q <= 1), 'Q must be in [0..1]');
  assert(all(diff(Q) > 0), 'Q must be strictly monotonically increasing');
  assert(all(diff(A) > -eps), 'A must be monotonically increasing');
  assert(all(diff(X) > -eps), 'X must be monotonically increasing');
  
  % Check extreme case where all values in X are outside range of A
  if all(X < A(1)+eps | X > A(end)-eps)
    F = nan(size(X));

  else
    % Extend Q and A at both ends
    Q = [0 Q(:)' 1];
    A = [-inf A(:)' inf];

    % Identify ranges where A is constant
    indAConst = find(diff(A) <= eps);
    k = find(diff([-inf indAConst]) > 1);
    firstIndAConst = indAConst(k);
    indAConst(k) = [];

    % Iterate over constant ranges
    noXCrit = false(size(firstIndAConst));
    for k = 1:length(firstIndAConst)
      ind = firstIndAConst(k);
      % Find largest values of X between current and preceding ranges
      xCrit = max(X(X < A(ind)-eps & X > A(ind-1)+eps));
      if isempty(xCrit)
        % No matching X found, corresponding (A,Q) may later be removed
        noXCrit(k) = true;
      else
        % Set (A,Q) to the linear interpolated value at xCrit
        if ind > 2
          Q(ind) = (xCrit-A(ind-1))/(A(ind)-A(ind-1))*(Q(ind)-Q(ind-1))+Q(ind-1);
        else
          Q(ind) = 0;
        end
        A(ind) = xCrit;
      end
    end

    % Remove remaining multiple points in A and correspondingly in Q, as well
    % as artificial extensions at both ends
    indRemove = [1 indAConst firstIndAConst(noXCrit) length(A)];
    A(indRemove) = [];
    Q(indRemove) = [];
    
    % If there are points in X that are by less than eps outside of the
    % range of A, adjust edges of A to allow interpolation at these X.
    xCrit = min(X(X < A(1) & X >= A(1)-eps));
    if ~isempty(xCrit)
      A(1) = xCrit;
    end
    xCrit = max(X(X > A(end) & X <= A(end)+eps));
    if ~isempty(xCrit)
      A(end) = xCrit;
    end

    % Invert A = A(Q) by linear interpolation
    F = interp1(A, Q, X, 'linear');
  end

  % Extrapolate at both ends with 0 and 1, resp.
  extrapolate = isnan(F);
  F(extrapolate & X < A(1)+eps) = 0;
  F(extrapolate & X > A(end)-eps) = 1;
end

