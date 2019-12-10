function [index, mzValueFound]=MSFindMZIndex(mzVector, mzValue)
% Return subscripts of the nearest spectrum value to mzValue in mzVector 
%   index=MSFindMZIndex(mzVector, mzValue)
%     mzVector: Input column-vector of mz-values (ascendent ordered)
%     mzValue:  Searched mzValue
%
%     index:    Index of nearest mz-value in mzVector to mzValue
% 
% For the result index it must hold:
% |mzVector(index)-mzValue|<=|mzVector(k)-mzValue| for all k

%some validation missing
if ~isnumeric(mzVector)||isempty(mzVector)
    error('mzVector must be a non-empty, numeric array');
end
if ~isnumeric(mzValue)||~isscalar(mzValue)||mzValue<0
    error('mzValue must be a non-negative scalar');
end
if(mzValue)<mzVector(1)
    index=1;
    return
end
if(mzValue)>mzVector(end)
    index=size(mzVector,2);
    return;
end
[~,index]=min(abs(mzVector-mzValue)); 
mzValueFound=mzVector(index);
