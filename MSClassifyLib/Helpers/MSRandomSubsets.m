function L = MSRandomSubsets (X, Nsub, F)
  % Generate random, disjoint subsets of specified sizes
  % L = MSRandomSubsets(N, Nsub): Generate a random vector L of length N
  %   with integer values>0 representing subsets of N items with sizes
  %   specified by the vector Nsub, In particular sum(L==j) == NSub(j) for 
  %   j=1..length(Nsub).
  % L = MSRandomSubsets(group, Nsub): As above, with additional group
  %   assignments given in the vector group, and the limitation that each 
  %   subset shall contain approximately the same proportion of elements in
  %   each group as the full set. Items without group are ignored.
  % L = MSRandomSubsets(group, Nsub, F): As above, with the vector F
  %   specifying the relative proportions of subset group sizes.
  %
  % TODO: The algorithm for selecting the random subsets is not very
  % sophisticated and is known to fail in case there are less items 
  % available in a particular class than required by the requested
  % class proportions. Should be replaced by a more robust algorithm.
  
  % Check input arguments
  narginchk(2,3);
  % First argument: either integer scalar N > 0 or group vector
  if isscalar(X)
    if ~(X > 0 && X == floor(X))
      error('Argument N must be a positive integer scalar');
    end
    N = X;
    % No group assignment to be used
    group = [];
    nGroups = 0;

  elseif isvector(X)
    % First argument is group vector
    group = X(:);
    N = sum(group > 0);
    if ~(N > 0 && all(group >= 0) && all(group == floor(group)))
      error(['Argument group must be a non-negative integer vector with '...
             'at least one non-zero element']);
    end
    nGroups = max(group);
    
  else
    error(['First argument must be a positive integer N, an integer vector ' ...
           'group or a scores matrix']);
  end
  
  % Second argument: vector of subset sizes
  if ~(~isempty(Nsub) && isvector(Nsub) && all(Nsub > 0 & floor(Nsub) == Nsub))
    error('Argument Nsub must be a vector of positive integers');
  elseif sum(Nsub) > N
    error('Total size of subsets exceeds number of items (=%d)', N);
  end
  Nsub = Nsub(:)';
  
  % Third argument: vector of relative output group proportions
  if nargin < 3 || isempty(F)
    % Not given, count group frequencies
    if isempty(group)
      F = 1;
    else
      M = group > 0;
      F = accumarray(group(M), 1);
    end
  elseif isempty(group)
    error('Argument F specified without group vector');
  else
    if ~(~isempty(F) && isvector(F) && all(F >= 0) && length(F) == nGroups)
      error(['Argument F must be a non-negative vector of length equal ' ...
             'the number of groups (=%d)'], nGroups);
    end
    F = F(:);
  end

  % Generate random ordering of numbers 1..N
  X = randperm(N);

  % Compute matrix of item counts for output subsets
  if nGroups <= 1
    % No group assignment and weights
    groupItems = {X};
    Nout = Nsub;

  else
    % Address items with group assignment
    I = find(group);
    X = I(X);
    % Split into items assigned to different groups
    groupItems = cell(nGroups,1);
    for j = 1:nGroups
      groupItems{j} = X(group(X) == j);
    end
    
    % Normalize F
    sF = sum(F);
    if sF > 0
      F = F(:)/sF;
    else
      F = ones(length(F),1)/sF;
    end
    Nout = zeros(nGroups,length(Nsub));
    for j = 1:nGroups-1
      Nout(j,:) = round(F(j)*Nsub);
    end
    Nout(end,:) = Nsub-sum(Nout,1);
  end
  
  % Set output label values
  if isempty(group)
    outLength = N;
  else
    outLength = length(group);
  end
  L = zeros(outLength,1);
  for j = 1:size(Nout,1)
    cumN = [0 cumsum(Nout(j,:))];
    for k = 1:length(Nsub)
      L(groupItems{j}((cumN(k)+1):cumN(k+1))) = k;
    end
  end
end

