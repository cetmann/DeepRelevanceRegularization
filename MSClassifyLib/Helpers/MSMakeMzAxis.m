function mzVector = MSMakeMzAxis (mzRange, D)
  % Generate m/z vector for given m/z range and bin count or width
  % mzVector = MSMakeMzAxis(mzRange, N): Generate m/z vector with N entries
  % mzVector = MSMakeMzAxis(mzRange, w): Generate m/z vector with bin width
  %   given by w. If w is scalar, bin width is constant. If w is a vector
  %   of length 2, bin width ranges from w(1) to w(2) (0 < w <= 1).
  %   mzRange: Positive vector of length 2
  
  % Check input arguments
  if ~(isnumeric(mzRange) && isvector(mzRange) && length(mzRange) == 2 && ...
       all(mzRange > 0) && all(diff(mzRange) > 0))
    error('mzRange must be a positive, increasing vector of length 2');
  end
  if ~(isnumeric(D) && isvector(D) && length(D) <= 2 && all(D > 0))
    D = [];  % D invalid: must be numeric, positive, scalar or 2-vector
  elseif length(D) == 2 && any(D > 1)
    D = [];  % D invalid: If 2-vector, it must be <= 1
  elseif length(D) == 1
    if D > 1
      if mod(D,1) > 0
        D = []; % D invalid: If scalar > 1, it must be integer
      else
        % Number of entries specified, convert to constant bin width
        D = diff(mzRange)/(D-1)*[1 1];
      end
    else
      % Constant bin width
      D = [D D];
    end
  end
  if isempty(D)
    error(['Second argument must be length of m/z vector, constant bin width, ' ...
           'or range of bin width values as [min max]']);
  end
  
  % Compute length of m/z vector
  N = round(2*(diff(mzRange)+D(2))/sum(D));
  % Compute bin width increment
  dD = diff(D)/(N-1);
  % Compute m/z vector
  n1 = 0:N-1;
  n2 = [0 0:N-2];
  mzVector = mzRange(1)+n1*D(1)+(n1.*n2)/2*dD;
end

