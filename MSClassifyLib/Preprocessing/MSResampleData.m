function resampledData = MSResampleData(data,mzVecOld,mzVecNew)
%resample data (matrix with spectra in rows) with associated mz vector
%mzVecOld onto new grid given by mzVecNew
%
%taking use of bioinformatic toolbox function msresample
%in case of down sampling (the recommended scenario!) before linear
%interpolation a low-pass filter is applied to reduce aliasing
%
%MSMaldiData has a method using MSResample
%
%in future this workaround might be changed such that a feature map can
%directly project onto mzVecNew

%check input arguments
if any([~isnumeric(data) ~isnumeric(mzVecOld) ~isnumeric(mzVecNew) ...
        ~ismatrix(data) ~isvector(mzVecOld) ~isvector(mzVecNew)])
    error(['data is expected to be a numeric matrix and mzVecOld ',...
        '/mzVecNew are expected to be numeric vectors']);
elseif length(mzVecOld) ~= size(data,2)
    error('data and mzVecOld dimension mismatch');
end

%choose windowType for the filter dependent on whether the 
%signal processing toolbox is available or not
verHelper = struct2cell(ver);
verHelper = squeeze(verHelper(1,1,:));
if any(strcmp(verHelper,'Signal Processing Toolbox'))
    windowType = 'kaiser';
else
%     warning('No signal processing toolbox found, choosing ''flattop'' window for the filter.');
    windowType = 'flattop';
end

%transposing mz vectors if necessary
mzVecOld = mzVecOld(:);
mzVecNew = mzVecNew(:);

[~,resampledData] = msresample(mzVecOld,data',mzVecNew,...
    'window',windowType,'cutoff',1);
%at the moment choosing nyquist frequency as cutoff
%(the maximum, so the filter has minimal effect)
resampledData = resampledData';
%during spline interpolation in msresample some negative values might occurr!!
resampledData(resampledData < 0) = 0;

end



