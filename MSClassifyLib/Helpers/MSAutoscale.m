function [data,dataMean,dataStd] = MSAutoscale(data,dataMean,dataStd)
if isempty(dataMean) && isempty(dataStd)
    dataMean = mean(data,1);
    dataStd = std(data,1);
end
data = bsxfun(@minus,data,dataMean);
data = bsxfun(@rdivide,data,dataStd);

% If one dimension is equals 0, dividing results in NaN entries.
% Set these values to 0.
data(isnan(data)) = 0;

end

