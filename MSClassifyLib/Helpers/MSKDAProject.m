% This helper function projects data to a (C-1)-dim. subspace of
% the corresponding reproducing kernel hilbert space, where C is the number
% of classes.
function features = MSKDAProject(data,projectionMatrix,kernel,...
    kernelParameter,trainingData)
features = (projectionMatrix' * ...
    MSKernelMatrix(kernel,trainingData',...
    data',kernelParameter))';

end