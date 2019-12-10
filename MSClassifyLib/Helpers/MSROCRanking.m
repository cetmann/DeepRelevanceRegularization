function [rankedList, rocMaxVals, rocValsTotal, pValsMax, pValsTotal, discrimClass] = MSROCRanking(msData,...
    labels, verboseFlag)
% Determine ROC values for each label vs the complement (of label>0)
% and in this way offer a ranking of the data
%
% Input:
% msData: The msData object can be both, an MSMaldiData object or an
%           MSFeatureData object
% labels: matching label data object (or a manually generated label vector)
% numOverallFeatures: (compare MSROCPicking) overall number of output
%           features
% verboseFlag: when set to true printing out the current calculation steps
%           (default: false)
%
% Output:
% All 3 output arguments are structs containing the 2 fields rocMax and
% rocMaxLabelBalance, in those fields they contain the following:
% rankedList: an index vector with the most significant features sorted to
%           to the top
% rocMaxVals: contains a vector with the corresponding highest ROC values 
%           among all labels (sorting matching to the one in rankedList)
% rocValsTotal: contains a matrix with the determined ROC values for all
%           (actually used!) labels (in the columns from left to right,
%           row sorting matching to the one in rankedList again)
% ---> All 3 output arguments are cropped to
%                 1:min(numOverall,msData.dataLength)
%

% Check input arguments
if ~isa(msData, 'MSData')
    error('msData must be an MSData object')
end
% Set defaults
if nargin < 3
    verboseFlag = false;
end
% Get label vector for specified classes
[L,nClasses,usedL] = MSLabelData.getLabelVector(labels, msData);
usedLGZero = usedL(usedL>0); %used labels~=0
numUsedLGZero = numel(usedLGZero);
% Must have at least two classes
if numUsedLGZero < 2
    error('At least two (non-empty) classes must be specified')
end
% In the special case of two classes iterating once (for the first class)
% is sufficient
if numUsedLGZero == 2
  labelIt = 1;
else
  labelIt = numUsedLGZero;
end

% Initialize array of one-vs-rest ROC values and p-values
nFeaturesAll = msData.dataLength;
rocValuesBasic = zeros(nFeaturesAll,labelIt+1);
pValuesBasic = zeros(nFeaturesAll,labelIt+1);
% Iterate over classes; determine ROC values and save them into
% rocValuesBasic
for k = 1:labelIt
    usedLk = usedLGZero(k);
    % Compute ROC values
    if verboseFlag
        fprintf(['Computing ROC values for class %d/%d (%d',...
            ' classes are non-empty)\n'], usedLk, nClasses, numUsedLGZero);
    end
    [r, p] = MSROC(msData.data(L == usedLk, :), ...
        msData.data(L > 0 & L~=usedLk, :),true);
    r=r';
    p=p';
    % Store top nFeatures values with corresponding indices and class
    % number into common list
    %[r,i] = sort(r, 'descend');
    rocValuesBasic(:,k) = r;
    pValuesBasic(:,k) = p;
end
% Add a column for the row-wise max values (can be used for rocMax-sorting)
[~, discrimClass] = max(abs(rocValuesBasic(:,1:labelIt)),[],2);
rocValuesBasic(:,labelIt+1) = maxAbsAux(rocValuesBasic, discrimClass);
pValuesBasic(:,labelIt+1) = maxAbsAux(pValuesBasic,discrimClass);
[~,rocMaxIndices] = sort(abs(rocValuesBasic(:,labelIt+1)),...
    'descend');

% rocMaxLabelBalance criterion
labelBestSortings = zeros(nFeaturesAll,labelIt);
labelBestSortingsROCs = zeros(nFeaturesAll,labelIt);
for k = 1:labelIt
    [labelBestSortingsROCs(:,k),labelBestSortings(:,k)] = sort(...
        abs(rocValuesBasic(:,k)),'descend');
end

% Row-wise vectorization
labelBestSortings = labelBestSortings';
labelBestSortingsROCs = labelBestSortingsROCs';
% write indices and values row-sorted into one matrix
rocMaxLabelBalIV = [labelBestSortings(:) labelBestSortingsROCs(:)];
for k = 1:nFeaturesAll
    rocMaxLabelBalIV((k-1)*labelIt+(1:labelIt),:) = ...
        sortrows(rocMaxLabelBalIV((k-1)*labelIt+(1:labelIt),:),-2);
end
[~, uniqueIndH] = unique(rocMaxLabelBalIV(:,1));
uniqueIndHSorted = sort(uniqueIndH);
rocMaxLabBalIndices = rocMaxLabelBalIV(uniqueIndHSorted);

% output arguments

rankedList.rocMax = rocMaxIndices;
rankedList.rocMaxLabelBalance = rocMaxLabBalIndices;
rocMaxVals = rocValuesBasic(:,labelIt+1);
rocValsTotal = rocValuesBasic(:,1:labelIt);

pValsMax = pValuesBasic(:,labelIt+1);
pValsTotal = pValuesBasic(:,1:labelIt);
end

function vector = maxAbsAux(A, indexes)
nrows = size(A,1);
vector = zeros(nrows,1);
for i=1:nrows
    vector(i)=A(i,indexes(i));
end
end




