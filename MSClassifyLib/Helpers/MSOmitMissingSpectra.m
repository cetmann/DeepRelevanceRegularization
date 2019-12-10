function [positions,spectraMatrix,idxD] = MSOmitMissingSpectra(positions,spectraMatrix)
%omitting missing spectra in raster matrix

% In cases where only very few spots on the raster matrix are missing,
% there is the option of simply omitting the corresponding spectra;
% this helper function transforms the position grid and the spectra matrix
% to fulfill exactly this purpose (deleting the respective rows in the
% spectraMatrix and "recounting" the indices in positions)

% the transformed positions and spectraMatrix are given as output along
% with idxD, a vector of the missing spectra indices,
%
% for direct generation of a MSMaldiData object via .sl-string this helper
% function is used in case spectra are missing in the raster matrix; in
% case later annotations are needed, the idxD vector is saved in the 
% property obj.metaInfo.omittedSpectraDuringImport
%

n = size(spectraMatrix,1);
idx = logical(positions);
idxD = find(~ismember(1:n , positions(idx)));

spectraMatrix(idxD,:) = [];
for k=1:length(idxD)
    idxA = idxD(k)-k+1;
    positions(positions>idxA) = positions(positions>idxA)-1;
end

fprintf('%g out of %g spectra (%.2f%%) removed!\n',numel(idxD),n,...
    100*numel(idxD)/n);

end

