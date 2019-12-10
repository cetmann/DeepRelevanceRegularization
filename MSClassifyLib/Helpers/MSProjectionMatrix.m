function projectionMatrix = MSProjectionMatrix(kernel,kernelParameter,...
  trainingData,classes,numFeatures,regularizationConstant)

numSamples = size(trainingData,1);
uniqueClasses = unique(classes)';
numClasses = length(uniqueClasses);

nClass = zeros(numClasses,1);
indicesClass = cell(numClasses,1);
mu_c = cell(numClasses,1);

K = MSKernelMatrix(kernel,trainingData(:,1:numFeatures)',[],kernelParameter);
mu = 1/numSamples * sum(K,2);
M = zeros(numSamples);
N = K*K';
for p = uniqueClasses
  if p~=0
    indicesClass{p} = (classes==p);
    nClass(p) = nnz(indicesClass{p});
    if nClass(p)>0
      Kc = MSKernelMatrix(kernel,trainingData(:,1:numFeatures)',...
        trainingData(indicesClass{p},1:numFeatures)',kernelParameter);
      mu_c{p} = (1/nClass(p)) * sum(Kc,2);
      M = M + nClass(p) * (mu_c{p} - mu)*(mu_c{p} - mu)';
      N = N - nClass(p) * (mu_c{p}*mu_c{p}');
    end
  end
end


[projectionMatrix,~] = ...
  eigs((N+regularizationConstant*eye(size(N)))\M,numClasses-1);

end
