function [Mnew,c] = MSHardThr(M,c)
%Helper: row-wise hard thresholding of input matrix

% M - matrix of spectra (row-wise)
% c - prctile value for the threshold (spectrum dependent)
% for each row of M all values below prctile(row,c) of that row are set 0
% (or the absolute values, to be more exact, negative peaks remain
% negative)
% 
% default value for c is set to 99.5
%

[n,~] = size(M);
Mnew = M;
if nargin == 1
    c = 99.5; %default thr parameter
end
for k = 1:n
    Mact = Mnew(k,:);
    Mact(abs(Mact)<prctile(abs(Mact),c)) = 0;
    Mnew(k,:) = Mact;
end

end


