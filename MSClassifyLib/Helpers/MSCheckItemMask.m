function [M, S] = MSCheckItemMask (D, M)
  % Check item mask or subset label and convert to logical matrix
  % M = MSCheckItemMask(D, M)
  %   D: Data object (MSData)
  %   M: Item mask, specified as an MSLabelData object, a logical vector
  %      or matrix, or as a numerical vector.
  %   If D is given as an empty array, size of M will not be checked and
  %   an empty M will not be converted to the default full item set.
  % [M, S] = MSCheckItemMask(D, M)
  %   Return names for item subsets in cell string array S. If M is not
  %   given as an MSLabelData object, default names '1', '2', ... are
  %   generated.

  if ~(isa(D, 'MSData') || isempty(D))
    error('Data argument must be an MSData object');
  end

  S = [];
  if ~isempty(M)
    % If M is non-empty, convert to logical matrix
    if islogical(M)
      % Item mask is given as a logical array
      if ~ismatrix(M)
        error('Logical item mask must be a vector or matrix');
      elseif isvector(M)
        % Logical vector, force to column vector
        M = M(:);
      end
    else
      if isa(M, 'MSLabelData')
        % Item mask is given as an MSLabelData object
        if M.dataLength ~= 1
          error('Item mask label data object must be single column');
        end
        S = M.labels;
        M = M.data;
      elseif isnumeric(M)
        if ~isvector(M)
          error('Numerical item mask must be a vector');
        end
        % Force to column vector
        M = M(:);
      else
        error('Item mask data type not supported');
      end
      % Convert numerical label vector to logical matrix
      M = accumarray([(1:length(M))' M+1], true, [], @any, false);
      M(:,1) = [];
    end
    % Check for empty item subsets
    if ~all(any(M,1))
      error('Item subsets must not be empty');
    end
  end
  % No label names? Generate default names
  if isempty(S)
    S = strtrim(cellstr(num2str((1:size(M,2))')))';
  end

  if ~isempty(D)
    if isempty(M)
      % By default, all items are used
      M = true(D.numItems,1);
    else
      % Check number of items
      if size(M,1) ~= D.numItems
        error('Item mask length must match number of data items');
      end
    end
  end
end
