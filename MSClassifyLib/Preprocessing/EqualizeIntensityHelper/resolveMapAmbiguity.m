function [X,Y] = resolveMapAmbiguity (X, Y, V, eps)
  % Resolve ambiguities in a map defined by (x,y) data points
  % [X,Y] = resolveMapAmbiguity(X,Y,V): Modify data points (X,Y) such that
  %   the resulting map is non-ambiguous and has the same effect on data
  %   values V.
  % [X,Y] = resolveMapAmbiguity(__, eps): Consider values with difference
  %   <= eps as equal (eps must be a non-negative scalar, default = 0).
  % X is assumed to be monotonically increasing.
  
  isVector = @(x) isnumeric(x) && isvector(x);
  if ~(isVector(X) && isVector(Y) && length(X) == length(Y))
    error('Arguments X and Y must be numeric vectors of same length');
  end
  if ~isnumeric(V) || isempty(V)
    error('Argument V must be a non-empty, numeric array');
  end
  if nargin < 4
    eps = 0;
  elseif ~(isnumeric(eps) && isscalar(eps) && eps >= 0)
    error('Argument eps must be a non-negative scalar');
  end
  
  % Find constant runs in X
  I = findConstantRuns(X, eps);
  
  if ~isempty(I)
    % Iterate over constant runs
    n = length(X);
    m = size(I,1);
    iDelete = false(n,1);
    for k = 1:m
      % Find closest data values in V falling into neighboring intervals
      i1 = I(k,1);
      vLo = max(V(V(:) < X(i1)-eps));
      if ~isempty(vLo) && i1 > 1 && vLo <= X(i1-1)
        vLo = [];
      end
      i2 = I(k,2);
      vHi = min(V(V(:) > X(i2)+eps));
      if ~isempty(vHi) && i2 < n && vHi <= X(i2+1)
        vHi = [];
      end

      if isempty(vLo) && isempty(vHi)
        % No data points in neighboring intervals, replace by center point
        Y(i1) = (Y(i1)+Y(i2))/2;
        iDelete(i1+1:i2) = true;
      elseif isempty(vLo) || vLo < (X(i1-1)+eps)
        % No data points in lower interval, remove lower map points
        iDelete(i1:i2-1) = true;
      elseif isempty(vHi) || vHi > (X(i2+1)-eps)
        % No data points in upper interval, remove upper map points
        iDelete(i1+1:i2) = true;
      else
        % Adjust lower and upper points by interpolation
        if i1 > 1
          m1 = (Y(i1)-Y(i1-1))/(X(i1)-X(i1-1));
          Y(i1) = Y(i1-1) + (vLo-X(i1-1))*m1;
        end
        if i2 < n
          m2 = (Y(i2+1)-Y(i2))/(X(i2)-X(i2+1));
          Y(i2) = Y(i2+1) + m2*(vHi-X(i2+1));
        end
        X(i1) = vLo;
        X(i2) = vHi;
        iDelete(i1+1:i2-1) = true;
      end

      % Remove obsolete map points
      X(iDelete) = [];
      Y(iDelete) = [];
    end
  end
end
