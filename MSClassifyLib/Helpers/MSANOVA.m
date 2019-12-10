function [fsNamesT,factorScoreT,highestComb,lowestComb,MEPamplitudes,...
    ANVAsignT,ANVAT,MEPYdata] = ...
    MSANOVA(factorM,scoreM,fsNames,continuousFactors,interactionFlag,displayFlag)
%calculate/display an ANOVA (analysis of variance) and a main effect plot
%
%Hint: ANOVA figures can be closed by typing "close all hidden"
%
%Interpretation hint: ANOVA tables: factors with Prob>F < 0.05 are
%   significant according to standard alpha = 0.05, with the input 
%   information for those factors it is unlikely they are simply a result 
%   of noise and random effects
%However, the 'significance' in ANOVA should not be overrated, usually the
%   standard interpretation of main effect plots is more reliable!
%Main effect plots: the general amplitude for the k-th subplot
%   corresponds to the effect a change of the k-th factor has on the
%   scores (according to the underlying linear model), on the y-axis the 
%   (in case of incomplete information estimated!) score means for fixed 
%   level of the k-the factor can be seen, if lower scores are "better" 
%   keep in mind that also the corresponding factor level that results in 
%   lower heights in these plots are estimated to be "better"
%
%Input: 
%factorM - matrix of integers describing different factor level
%   combinations, the number of columns is the number of factors which are
%   varied, the number of rows is the number of tested combinations (not 
%   all possible combinations have to be included);
%   the integers present the factor levels, identical integers
%   correspond to the same value;
%   MSANOVAgenFactorMatrix.m offers help in generating a (complete) version of
%   factorM
%scoresM - column vector / matrix of scores for the factor level 
%   combinations, the number of rows has to be the same as in factorM;
%   the number of columns can be greater than one, then the analysis is
%   repeated for each column of scoresM
%Optional input parameters:
%fsNames - optional: if given as input (fsNames = [] is treated like no 
%   input), it is expected to be a cell array
%   of strings, the number of columns is expected to be the same as the sum
%   of columns in factorM and scoreM, the number of rows is expected to be
%   as high as the maximum number of factor levels for any factor plus 1;
%   the first row contains names for the factors / scores,
%   the second and following row names for corresponding factor levels
%   (in a monotone matching order to the integers in factorsM)
%   factorM will be transformed (monotonically increasing) such that for
%   each column only the integers 1:m exist where m is the number of
%   factor levels (can vary from column to column) - e.g. psNames(3,4)
%   contains the name of the 2nd factor level for the fourth factor
%   unused entries have to be assigned the empty string '' for 
%   flexibility reasons in factorM; (compare fsNamesDummy)
%continuousFactors - optional: factors can be categorical (like
%   nationality, gender, ...) or continuous (like time, temperature, ...);
%   specifying which factors are continuous in the vector continuousFactors
%   (e.g. continuousFactors = [2 5] if factors #2 and #5 are "continuous")
%   improves the approximation to reality of the underlying statistical 
%   model; by default all factors are assumed to be categorical, i.e.
%   continuousFactors = [];
%interactionFlag - binary flag (default is 0); if set to 1 the ANOVA model
%   is changed to the 'interaction' model instead of the default 'linear'
%   model
%displayFlag - showing graphs if true (default is true)
%
%Output: 
%fsNamesT - fsNames in a table form,
%factorScoreT - scoreM in combination with factorM in table form
%highestComb - for each score the combination with highest scores
%               (according to the linear model) in table form
%lowestComb - for each score the combination with lowest scores
%               (according to the linear model) in table form
%MEPamplitudes - amplitudes of the main effects in table form
%ANVAsignT - summary of significant factors according to the ANOVA
%           (alpha<0.05)
%ANVAT - (displayed) ANOVA tables for each score in a cell array
%MEPYdata - (displayed) MEP y data in a cell array


%define important dimensions of input variables
numFactors = size(factorM,2);
numCombinations = size(factorM,1);
numScoreTypes = size(scoreM,2);

%check numeric type of input variables factorM and scoreM
if ~isnumeric(factorM) || ~isnumeric(scoreM)
    error('factorM and scoreM have to be numeric!');
end

%check input displayFlag, interactionFlag, continuousFactors
if nargin < 6
  displayFlag = true;
end
if displayFlag
  displayFlagS = 'on';
else
  displayFlagS = 'off';
end
if nargin < 5
    interactionFlag = 0;
end
if nargin < 4
    continuousFactors = [];
end
if ~isnumeric(continuousFactors) || ...
        (~isempty(continuousFactors) && (~isvector(continuousFactors) ...
         || norm(abs(round(continuousFactors))-continuousFactors)>0 ...
         || min(continuousFactors)<1 || max(continuousFactors)>numFactors))
    error(['continuousFactors has to be a vector of integers in the',...
        ' range of 1:numFactors!']);
end

%check if factorM contains only integers
if round(factorM) ~= factorM
    error('factorM does contain non-integers!');
    %the fact itself is unproblematic but the danger of confusion of some
    %factorM and scoreM columns is important enough to include this error
end

%transform integers in factorM to 1:m-form
numFactorLvls = zeros(1,numFactors);
factorMnew = zeros(size(factorM));
for k = 1:numFactors
    factorLvlsK = unique(factorM(:,k));
    numFactorLvls(k) = numel(factorLvlsK);
    for l = 1:numFactorLvls(k)
        factorMnew( factorM(:,k)==factorLvlsK(l),k )=l;
    end
end

%consistency check of factor levels with fsNames
if any(sum(cellfun(@(x)~isempty(x),fsNames(2:end,1:numFactors)),1)-...
    numFactorLvls)
    error(['The number of actually occurring factor levels is',...
        'inconsistent with fsNames!']);
end
factorM = uint32(factorMnew);

%check dimension match of scoreM and factorM - the number of rows has to be
%identical
if size(scoreM,1) ~= numCombinations
    error('scoreM has to have the same number of rows as factorM');
end

%check cellstr type of optional input variable fsNames and dimensions
if nargin >= 3 && ~isempty(fsNames) 
    if ~(iscellstr(fsNames) || iscell(fsNames)) || ...
        size(fsNames,2)~=numFactors+numScoreTypes...
            || size(fsNames,1) < max(numFactorLvls)+1
        error(['fsNames has to be a cell array of strings with',... 
            ' size(factorM,2)+size(scoreM,2) columns and (at least)',...
            ' max(numFactorLvls)+1 rows']);
    end
else
    fsNames = genfsNamesDummy(numFactorLvls,numScoreTypes);
end
if iscell(fsNames{1,2})
  fsNamesH = fsNames;
  ccF = cellfun(@(x)numel(x),fsNames(:,2));
  fsNames = cell( size(fsNamesH,1), 1+max(ccF) );
  fsNames(:,1) = fsNamesH(:,1);
  for kk = 1: size(fsNamesH,1)
    fsNames{:,2:(1+ccF(kk))} = fsNamesH{kk,2};
  end
end


fsNamesT = array2table(fsNames(2:end,:),'VariableNames',fsNames(1,:),...
    'RowNames',cellstr([char(repmat('factorLvl',max(numFactorLvls),1)) ...
    num2str((1:max(numFactorLvls))')]));


%--------------------------------------------------------------------------
%---------------------- main part -----------------------------------------

% showing ANOVA tables and main effect plots
ANVA = cell(1,numScoreTypes);
ANVAT = cell(1,numScoreTypes);
MEPfigH = cell(1,numScoreTypes);
MEPaxH = cell(1,numScoreTypes);
MEPChildren = cell(1,numScoreTypes);
MEPYdata = cell(numFactors,numScoreTypes);

for k = 1:numScoreTypes
    if interactionFlag
        [ANVA{k},ANVAT{k}] = anovan(scoreM(:,k),factorM,'varnames',...
            fsNames(1,1:numFactors),'continuous',continuousFactors,...
            'model','interaction','display',displayFlagS);
    else
        [ANVA{k},ANVAT{k}] = anovan(scoreM(:,k),factorM,'varnames',...
            fsNames(1,1:numFactors),'continuous',continuousFactors,...
            'model','linear','display',displayFlagS);
    end
    titleA = findall(0,'Type','figure');
    titleA(1).Name = ['ANOVA for ' fsNames{1,numFactors+k}];
    
    if displayFlag
      figure('Name',['MEP for ' fsNames{1,numFactors+k}]);
      [MEPfigH{k},MEPaxH{k}] = maineffectsplot(scoreM(:,k),...
          factorM,'varnames',fsNames(1,1:numFactors));
      MEPChildren{k} = get(MEPaxH{k},'children');
      for l = 1:numFactors
          MEPYdataHelper = MEPChildren{k};
          MEPYdata{l,k} = MEPYdataHelper{l}.YData;
          %setting xTickLabels manually
          xTickHelper1 = MEPaxH{k}(l);
          xTickHelper2 = xTickHelper1.get.XTickLabel;
          xTickHelper1.XTickLabel = fsNames(1+...
              str2double(xTickHelper2),l);
          MEPaxH{k}(l).XTickLabelRotation = 70;
      end
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%create table of input data 
%(in combination with fsNames it contains the information in a neat form)
factorScoreT = array2table([double(factorM) double(scoreM)],...
    'VariableNames',fsNames(1,:));

%determine special combinations with highest and lowest score (from the
%estimated distribution, there does not have to exist a
%measurement/experiment for the exact combination!)
[maxCombVals,maxCombHelp] = cellfun(@max,MEPYdata);
highestComb = cell(size(maxCombHelp));
for l = 1:numFactors
    for k = 1:numScoreTypes
        highestComb{l,k} = fsNames{1+maxCombHelp(l,k),l};
    end
end
highestComb = array2table(highestComb,'VariableNames',...
    fsNames(1,(numFactors+1):end),'RowNames',...
    fsNames(1,1:numFactors));
[minCombVals,minCombHelp] = cellfun(@min,MEPYdata);
lowestComb = cell(size(minCombHelp));
for l = 1:numFactors
    for k = 1:numScoreTypes
        lowestComb{l,k} = fsNames{1+minCombHelp(l,k),l};
    end
end
lowestComb = array2table(lowestComb,'VariableNames',...
    fsNames(1,(numFactors+1):end),'RowNames',...
    fsNames(1,1:numFactors));

%calculate amplitudes of MEP
MEPamplitudes = maxCombVals - minCombVals;
MEPamplitudes = array2table(MEPamplitudes,'VariableNames',...
    fsNames(1,(numFactors+1):end),'RowNames',...
    fsNames(1,1:numFactors));

%transform ANVA (contains Prob>F values) to table
% ANVAT = array2table(cell2mat(ANVA),'VariableNames',...
%     fsNames(1,(numFactors+1):end),'RowNames',...
%     fsNames(1,1:numFactors));

rowNamesInteractions = ANVAT{1};
rowNamesInteractions = rowNamesInteractions(2:(end-2),1);

%extract significant factors
ANVAsignhelp = cell2mat(ANVA)<0.05; %standard level of significance
ANVAsign = cell(size(ANVAsignhelp));
ANVAsign(:) =               {'     -     '};
ANVAsign(ANVAsignhelp) =    {'significant'};
ANVAsignT = array2table(ANVAsign,'VariableNames',...
    fsNames(1,(numFactors+1):end),'RowNames',...
    rowNamesInteractions);

end


function fsNamesDummy = genfsNamesDummy(numFactorLvls,numScoreTypes)
%generating dummy of fsNames; useful if there is no fsNames input in
%MSANOVA

%input: numFactorLvls - 1xn vector (for n factors = n columns in factorM),
%           the k-th entry corresponds to the number of different factor 
%           levels for the k-th factor
%       numScoreTypes - number of different score types = number of columns
%           in fsNames
%

fsNamesDummy = cell(max(numFactorLvls)+1,...
    length(numFactorLvls)+numScoreTypes);
fsNamesDummy(:) = {''};
for k = 1:length(numFactorLvls)
    fsNamesDummy{1,k} = ['factor' num2str(k)];
    for l=1:numFactorLvls(k)
        fsNamesDummy{1+l,k} = ['factor' num2str(k) 'Lvl' num2str(l)];
    end
end
for k = 1:numScoreTypes
    fsNamesDummy{1,length(numFactorLvls)+k} = ['scoreType' num2str(k)];
end

end

