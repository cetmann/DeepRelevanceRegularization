function [Mnew,c] = MSSoftThr(M,c)
%Helper: row-wise soft thresholding of input matrix
%
% M - matrix of spectra (row-wise)
% c - prctile value for the threshold (spectrum dependent)
%  (determines threshold value)
% 
% default value for c is set to 99.5
%

[n,~] = size(M);
Mnew = M;
if nargin == 1
    c = 99.5; %default thr parameter
end
for k = 1:n
    Mcurr = Mnew(k,:);
    currRowThr = prctile(abs(Mcurr),c);
    Mnew(k,:) = wthresh(Mcurr,'s',currRowThr);
end

end
