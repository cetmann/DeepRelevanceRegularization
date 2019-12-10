function cellNVPairsA = struct2NameValuePair( structA )
% Transforms a struct to a name value pair cell, used e.g. MSSVMClassifier
  f = fields( structA );
  n = numel(f);
  if n == 0
    cellNVPairsA = [];
  else
  
    cellNVPairsA = cell( 1, 2*n );
    for i = 1:n
      idx1 = 2*(i-1)+1;
      idx2 = idx1+1;
      cellNVPairsA{ idx1 } = f{i};
      cellNVPairsA{ idx2 } = structA.(f{i});
    end
    
  end
end