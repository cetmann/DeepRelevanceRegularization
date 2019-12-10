function factorM = MSANOVAgenFactorMatrix(numFactorLvls)
%generate a factorM matrix that can be used directly or postprocessed for 
%MSANOVA
%
%The type of matrix generated does not contain repetitions, all possible
%combinations of factor levels occur once
%
%Input (numFactorLvls) is a vector containing the number of factor levels, 
%e.g. if 3 factors are varied, where the first and third have 2 
%and the second 3 levels, 
%the input of choice would be numFactorLvls = [2 3 2];
%less than 2 or more than 9 factor levels for single factors are not 
%supported

if ~isnumeric(numFactorLvls) || ...
        norm(double(uint32(numFactorLvls))-double(numFactorLvls))>0 ...
        || any(numFactorLvls < 2) || any(numFactorLvls > 9)
    error('numFactorLvls can only contain integers > 1 and < 10');
end
m = max(numFactorLvls);
l = length(numFactorLvls);
factorM = dec2base(0:(m^l-1),m);
factorM = cellfun(@(row)str2double(regexp(num2str(row),'\d','match')),...
    cellstr(factorM),'UniformOutput',0);
factorM = 1+cell2mat(factorM);
for k = 1:l
    factorM(factorM(:,k)>numFactorLvls(k),:) = [];
end
factorM = uint8(factorM);

end

