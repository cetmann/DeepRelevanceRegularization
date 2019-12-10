function T = MSContingencyTable (L1, L2, S)
  % Generate contingency table from label data
  % T = MSContingencyTable(L1, L2): Generate a table T whith rows
  %   and columns corresponding to label values in L1 and L2, resp. Table
  %   entries count the number of items with the corresponding combination
  %   of label values.
  % T = MSContingencyTable(L): Generate a single column table listing the
  %   frequencies of the individual label values in L.
  % T = MSContingencyTable(L, varName): Use varName as the name for the
  %   frequency table, default is "Count".
  % T = MSContingencyTable(__, S): Perform 'unique' aggregation to label
  %   objects before computing the contingency table.
  %   Note: Use MSContingencyTable(L, [], S) to omit second label argument
  % L, L1, L2, S must be MSLabelData objects with a single data column
  
  % Check input arguments
  narginchk(1,3);
  % Default variable name for single column output
  varNames = {'Count'};
  % Check optional argument
  if nargin < 2
    L2 = [];
  elseif ischar(L2)
    varNames = {matlab.lang.makeValidName(L2)};
    L2 = [];
  end
  if nargin < 3
    S = [];
  end
  % Check label objects
  if ~(isa(L1, 'MSLabelData') && L1.dataLength == 1) || ...
     ~(isempty(L2) || (isa(L2, 'MSLabelData') && L2.dataLength == 1)) || ...
     ~(isempty(S) || (isa(S, 'MSLabelData') && S.dataLength == 1))
    error('Label input arguments must be MSLabelData objects with one data column');
  end

  % Perform aggregation
  if ~isempty(S)
    L1 = AggregateLabel(L1, S);
    if ~isempty(L2)
      L2 = AggregateLabel(L2, S);
    end
  end
  
  % Prepare table computation
  M = L1.data > 0;
  rowNames = matlab.lang.makeValidName(L1.labels);
  
  if isempty(L2)
    % Frequency table for L1
    T = accumarray(L1.data(M), 1, [L1.numLabels 1]);
  else
    % Contingency table for (L1,L2)
    M = M & L2.data > 0;
    T = accumarray([L1.data(M), L2.data(M)], 1, [L1.numLabels, L2.numLabels]);
    varNames = matlab.lang.makeValidName(L2.labels);
  end
  T = array2table(T, 'RowNames', rowNames, 'VariableNames', varNames);
end

function L = AggregateLabel (L, S)
  % Aggregate labels in L with respect to S and add a 'MIXED_' label in
  % case unique aggregation fales
  L = L.aggregate(S, 'unique');
  if any(isnan(L.data(:)))
    L.setLabels([L.labels; 'MIXED_']);
    L.data(isnan(L.data(:))) = L.numLabels;
  end
end