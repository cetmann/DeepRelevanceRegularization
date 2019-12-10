function MSSmooth (msData)
% Perform smooting of spectra with (1,2,1) weighted average

% Check argument
if ~isa(msData, 'MSData')
  error('Argument must be an MSData object')
end

% Compute weighted average spectra
fill = zeros(msData.numItems,1);
msData.data = ([fill,msData.data(:,1:end-1)]+[msData.data(:,2:end),fill] ...
                  + 2*msData.data)/4;

end

