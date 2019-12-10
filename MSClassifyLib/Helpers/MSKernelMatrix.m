% Function that efficiently computes kernel matrices.
% If Y is empty, the gram matrix (with respect to the specific kernel
% function) is calculated. Otherwise, the (i,j)-th entry of the kernel
% matrix is k(x_i,y_j).
function K = MSKernelMatrix(kernel,X,Y,kernelParameter)
if strcmp(kernel,'rbf')
  if ~isempty(Y)
    X = double(X);
    Y = double(Y);
    [~,numXSamples] = size(X);
    numYSamples = size(Y,2);
    K = (ones(numYSamples,1)*sum(X.^2,1))';
    K = K + (ones(numXSamples,1)*sum(Y.^2,1));
    K = K - 2*(X'*Y);
    K = exp((-kernelParameter)*K);
  else
    X = double(X);
    [~,numXSamples] = size(X);
    K = (ones(numXSamples,1)*sum(X.^2,1))';
    K = K + ones(numXSamples,1)*sum(X.^2,1);
    K = K - 2*(X'*X);
    K = exp((-kernelParameter)*K);
  end
end
if strcmp(kernel,'linear')
  if ~isempty(Y)
    Y = double(Y);
    X = double(X);
    K = X'*Y;
  else
    X = double(X);
    K = X'*X;
  end
end
if strcmp(kernel,'quadratic')
  if ~isempty(Y)
    Y = double(Y);
    X = double(X);
    K = (X'*Y).^2;
  else
    X = double(X);
    K = (X'*X).^2;
  end
end
if strcmp(kernel,'sigmoid')
  if ~isempty(Y)
    Y = double(Y);
    X = double(X);
    K = X'*Y;
    K = tanh(K + kernelParameter * ones(size(K)));
  else
    X = double(X);
    K = X'*X;
    K = tanh(K + kernelParameter * ones(size(K)));
  end
end

end