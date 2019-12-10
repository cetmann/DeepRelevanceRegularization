function indA2B = MSResorter(A,B)
% If there is an original sorting (1:n)' and a sorting A and a sorting B
% (vectors of length n containing a permutation of indices 1:n) 
% this function has a function handle as an output that resorts the rows
% of an input (n,:) input-matrix from A to B (current sorting of matrix is
% according to A, output sorting is in B)
%
% is used in MSLinearBasisMap for resorting

[~,sortIndA] = sort(A);
indA2Origin = @(x)x(sortIndA,:);
indOrigin2B = @(x,B)x(B,:);
indA2B = @(x)indOrigin2B(indA2Origin(x),B);

end


