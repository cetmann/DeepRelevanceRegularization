function MSAlign (maldiData, nRefPeaks, minMzDist)
% Perform alignment of spectral peaks on all spectra in dataset
% Check argument
narginchk(1,3)
if ~isa(maldiData, 'MSMaldiData')
  error('Argument must be an MSMaldiData object')
end

if nargin<3 || isempty(minMzDist); 
    minMzDist = 1.5; % Dalton, default value
elseif ~isscalar(minMzDist) || ~isnumeric(minMzDist) || minMzDist<0
        error('minMzDist must be a real non-negative number')
end

if nargin<2 || isempty(nRefPeaks)
    nRefPeaks = 10; %default value
else
    if ~isscalar(nRefPeaks) || ~isnumeric(nRefPeaks) || nRefPeaks<1
        error('minMzDist must be a positive integer')
    end
    nRefPeaks=floor(nRefPeaks);
end



% Compute mean spectrum
meanSpec = mean(maldiData.data, 1);

% Find reference peaks in mean spectrum
refPeaks = MSFindLargestPeaks(meanSpec, maldiData.mzVector, nRefPeaks, minMzDist);

% Compute aligned spectra
msalignParams = {'MaxShift',[-1 1]; 'WidthOfPulses', maldiData.mzResolution}';
maldiData.data = msalign(maldiData.mzVector', maldiData.data', ...
                         maldiData.mzVector(refPeaks), msalignParams{:})';
maldiData.data(isnan(maldiData.data)) = 0;

end



