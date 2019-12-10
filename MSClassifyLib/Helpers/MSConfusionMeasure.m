function Y = MSConfusionMeasure (X, measure)
  % Compute classification performance measures from a confusion matrix
  % Y = MSConfusionMeasure(X): Compute several performance measures from
  %   an n-by-n input matrix X (n >= 2) and return results in a struct Y.
  %
  % TODO: More detailed description
  
  % Check input arguments
  narginchk(1,2);
  if ~(isnumeric(X) && size(X,1) == size(X,2) && size(X,1) >= 2 && ...
       all(X(:) >= 0))
    error('Argument must be a numeric, square matrix of non-negative integers');
  end
  
  % Size of input matrix. ND-arrays (N>2) are treated as arrays of square
  % matrices.
  sx = size(X);
  numClasses = sx(1);
  
  % True and predicted per class totals
  predTotal = sum(X,1);
  trueTotal = reshape(sum(X,2), size(predTotal));
  % Overall total count
  total = sum(trueTotal,2);
  
  % Correct per class prediction counts
  correct = zeros([1, sx([1, 3:end])]);
  for k = 1:numClasses
    correct(1,k,:) = X(k,k,:);
  end
  
  % Per class performance measures
  sensitivity = correct./trueTotal;
  specificity = 1-(predTotal-correct) ./ ...
                  (repmat(total, [1, numClasses])-trueTotal);
  
  % Multi-class performance measures
  accuracy = sum(correct,2)./total;
  kappa = (accuracy*numClasses-1)/(numClasses-1);
  
  % excluding NaN cases to generate avg. sensitivy / specificity 
  nansSensitivity = isnan(sensitivity);
  numNonNansSe = sum(~nansSensitivity,2);
  nansSpecificity = isnan(specificity);
  numNonNansSp = sum(~nansSpecificity,2);
  sensitivityHelper = sensitivity;
  specificityHelper = specificity;
  sensitivityHelper(nansSensitivity) = 0;
  specificityHelper(nansSpecificity) = 0;
  
  avgSensitivity = sum(sensitivityHelper,2)./numNonNansSe;
  avgSpecificity = sum(specificityHelper,2)./numNonNansSp;
  
  % Create result struct, remove singletons in first dimension
  Y = struct('accuracy',       shiftdim(accuracy, 1), ...
             'kappa',          shiftdim(kappa, 1), ...
             'sensitivity',    shiftdim(sensitivity, 1), ...
             'specificity',    shiftdim(specificity, 1), ...
             'avgSensitivity', shiftdim(avgSensitivity, 1), ...
             'avgSpecificity', shiftdim(avgSpecificity, 1));
  
  if nargin >= 2
    Y = Y.(measure);
  end
end

