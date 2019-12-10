function A = spatialKernelToMatrix (kernel, positions)
  % Check arguments
  if ~(isnumeric(kernel) && ismatrix(kernel) &&~isempty(kernel))
    error('Kernel argument must be a non-empty, numeric matrix');
  elseif abs(sum(kernel(:))) < 1e-7
    error('Sum of kernel elements is zero, not supported');
  end
  if ~isa(positions, 'MSPositionGrid')
    error('Positions argument must be valid MSPositionGrid object');
  end
  
  % Compute rows and columns for each grid item and kernel element
  m = positions.numItems;
  n = numel(kernel);
  [itemRow, itemCol] = ind2sub(positions.gridSize, positions.reverseIndex);
  [kernelRow, kernelCol] = ind2sub(size(kernel), 1:n);
  kernelCenter = floor((size(kernel)-1)/2)+1;
  rowMatrix = repmat(itemRow, [1 n])+repmat(kernelRow, [m 1])-kernelCenter(1);
  colMatrix = repmat(itemCol, [1 n])+repmat(kernelCol, [m 1])-kernelCenter(2);
  
  % Compute item indices for each grid item and kernel element
  rowMatrix(rowMatrix < 1 | rowMatrix > positions.gridSize(1)) = nan;
  colMatrix(colMatrix < 1 | colMatrix > positions.gridSize(2)) = nan;
  indMatrix = sub2ind(positions.gridSize, rowMatrix, colMatrix);
  indValid = ~isnan(indMatrix);
  indMatrix(indValid) = positions.indexGrid(indMatrix(indValid));
  indValid = indMatrix > 0;
  
  % Generate convolution matrix entries
  cnvMatrix = repmat(kernel(:)', [m 1]) .* indValid;
  % Rescale entries to account for kernel positions outside boundary
  kernelSum = sum(kernel(:));
  cnvMatrix = cnvMatrix .* repmat(kernelSum./sum(cnvMatrix,2), [1 n]);
  
  % Generate sparse convolution matrix
  Arow = repmat((1:m)', [1 n]);
  A = sparse(Arow(indValid), indMatrix(indValid), cnvMatrix(indValid), m, m);
end

