function peaks = MSFindLargestPeaks(spc, mzVector, nLargest, minDist)
% Return subscripts of the largest peaks in a spectrum
%   peaks = findLargestPeaks(spc, mzVector, nLargest, minDist)
%     spc:      Input Spectrum
%     mzVector: Vector of mz-values, same length as spc
%     nLargest: Number of largest peaks to identify
%     minDist:  Minimum distance between two large peaks (in mz-units, default = 1)
%
%     peaks:    Sorted subscripts of spectrum peaks
% 
% Identify peaks that dominate surrounding spectrum intensities.
% For a subscript K to be included in peaks, the following must hold:
%   spc(K) > spc(I) for all I with |mzVector(I)-mzVector(K)| < minDist

if nargin < 4
  minDist = 1;
end

[~,K] = sort(spc(:),'descend');

peaks = [];

if nLargest >= 1 && length(K) >= 1
  peaks = K(1);

  if nLargest > 1
    for i = 2:length(K)
      if min(abs(mzVector(K(i))-mzVector(K(1:i-1)))) >= minDist
        peaks(end+1) = K(i);
        if length(peaks) >= nLargest
          break
        end
      end
    end
  end
end

end
