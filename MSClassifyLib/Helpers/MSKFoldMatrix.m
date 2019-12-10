function kFold = MSKFoldMatrix(labelVector)
    % This helper function lets you represent a label vector as a
    % kFold-matrix, where 1 stands for a member of the training set and a
    % 2 stands for a member of the test set (of a certain CV-fold).

    % First, one-hot encode the labelVector. Leaves 0 as 0.
    oneHotEncodedLabels = bsxfun(...
        @eq,...
        labelVector(:),...
        1:max(labelVector));
    nCols = size(oneHotEncodedLabels,2);
    
    % Add 1 where there is a label.
    onesMatrix = repmat(sum(oneHotEncodedLabels,2),1,nCols);
    kFold = oneHotEncodedLabels + onesMatrix;
end
    
    