function patternCellOut = reducePatternCell(patternCell, minim, maxim)
% Helper function which selects randomly a few peptides of the patternCell
% protein to return a new 'expressed' protein by those chosen peptides.
%   INPUT
%       patternCell: MSPaternCell object representing the targeted proteins
%       minim: minimum number of peptides chosen
%       maxim: maximum number of peptides chosen
if ~isa(patternCell,'MSPatternCell')
    error('The first argument must be of type MSPatternCell')
end
if ~isPositiveInteger(minim)||~isPositiveInteger(maxim)
    error('Second and third arguments must be positive integers')
end
patternCellOut = cell(1,patternCell.nPatterns);
for i=1:patternCell.nPatterns
    p = patternCell.patterns{i};
    rows = size(p,1);    
    n = randi([minim, maxim]);
    n = min([rows,n]);
    patternCellOut{i} = p(datasample(1:rows,n,'Replace',false),:);
end
patternCellOut = MSPatternCell(patternCellOut,patternCell.mzVector,patternCell.proteinNames);
end
function bool = isPositiveInteger(number)
    bool = isscalar(number)&& isnumeric(number) && number>0 && floor(number)==number;
end