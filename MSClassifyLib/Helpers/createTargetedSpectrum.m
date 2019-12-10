function peptidPattern = createTargetedSpectrum(mzVector, filename, taskName, numberSheets)
% INPUT
%  mzVector: vector of mz-values for which the pattern will be generated
%  filename: path+name of the excell tables containing information of
%  peptides
%  taskName: (optional) Identifier of the dataset for which the biomarkers are relevant   
%  (filename for OvsB: '/home/yovany/Digipath/TargetedNMF/Biomarker - MSI
%  verified and SOP - BR OV.xlsx')
%  (filename for ADSQ: '/home/yovany/Digipath/TargetedNMF/Biomarker -
%  ADvsSQ.xlsx')
narginchk(2,4)
if nargin < 3
    taskName = [];
end
if nargin < 4
    numberSheets = 0;
end


% Create spectra for a biomarker given by a list of peptides:

%% specify spectrometry specific parameters and mz-axis
fwhh =  2932/18000; % Taken from MIL-Homepage (instrumentation) R = M / dM
nMz = length(mzVector);
%% Specific chemical and table parameters
% specify chemical parameters
% protonMass = 1.008;
% oxidMass = 15.9949;
% methioninMass = 131.0404;

% specific parameters of tables
%indexStart = 5;
indexStart = 12;
modificationIndex = 4;
sequenceIndex = 8;
%% Load list of peptides and create patterns

% for getting name of sheets
[~,nameProteins]=xlsfinfo(filename);
if numberSheets == 0
    numberSheets = length(nameProteins);
end
nameProteins = nameProteins(1:numberSheets);
%numberSheets = length(nameProteins);
peptidQ=cell(1,numberSheets);
for sheet = 1: numberSheets
    % read excel table
    [num,txt,raw] = xlsread(filename, sheet);
    % determine last index (in the last sheet there are additional rows)
    for i = indexStart:size(raw,1)
        if raw{i,1} ~= 1
           i = i-1;
           break; 
        end
    end
    lastIndex = i;

    % extract relevant columns from table
    rawCell = cell(lastIndex - indexStart + 1,2);
    rawCell(:,1) = raw(indexStart:lastIndex, modificationIndex);
    rawCell(:,2) = raw(indexStart:lastIndex, sequenceIndex);

    % Create isotopic patterns separately
    lenList = size(rawCell,1);
    combindedPattern = zeros(1, size(mzVector, 2));
    Q=zeros(lenList,nMz);
    insideInterval = false(lenList,1);
    for j = 1: lenList
%         queryStr = rawCell{j,1};
%         rawSequence = rawCell{j,2};
%         while strcmp(rawSequence(end),' ')
%             rawSequence = rawSequence(1:end-1);
%         end
%         sequence = rawSequence(4: end-3);
%         % Acetly modifications
%         if strfind(queryStr, 'Acetyl')  
%            [MD, Info, DF] = isotopicdist(sequence, ...
%                 'Resolution', fwhh, 'NTerminal', 'acetyl');
%         else  % without acetyl modification
%            [MD, Info, DF] = isotopicdist(sequence, ...
%                'Resolution', fwhh);
%         end
%         % Oxidation
%         if strfind(queryStr, 'Oxidation')
%             numOxidations = str2double(queryStr(...
%                 strfind(queryStr, 'Oxidation') - 1));
%             DF(:,1) = DF(:,1) + numOxidations * oxidMass;
%         end
%         % Met-loss
%         if strfind(queryStr, 'Met-loss')
%             numMetloss = str2double(queryStr(...
%                 strfind(queryStr, 'Met-loss') - 1));
%             DF(:,1) = DF(:,1) - numMetloss * methioninMass;
%         end
%         % add H+ Ion (translate by 1.008 Da to the right)
%         DF(:,1) = DF(:,1) + protonMass;
        queryStr = rawCell{j,1};
        rawSequence = strtrim(rawCell{j,2});
        sequence = rawSequence(4: end-3);
        % Met-loss (remove aminoacid M from beginning of sequence)
        if strfind(queryStr, 'Met-loss')
            sequence = sequence(2:end);
        end
    
        % get structure formula from amino-acid sequence
        [~, Info, ~] = isotopicdist(sequence, 'Resolution', fwhh);
        
        % Oxidation (add oxid atom(s) to the structure formula)
        if strfind(queryStr, 'Oxidation')
            numOxidations = str2double(queryStr(...
                strfind(queryStr, 'Oxidation') - 1));
            Info.Formula.O = Info.Formula.O + numOxidations;
        end
        
        % add H+ Ion to the structure formula 
        Info.Formula.H = Info.Formula.H +1;
        
        % recompute isotopes (and check for acetly modifications)
        if strfind(queryStr, 'Acetyl')  
            [MD, Info, DF] = isotopicdist(Info.Formula, ...
                    'Resolution', fwhh, 'NTerminal', 'acetyl');
        else  % without acetyl modification
        	[MD, Info, DF] = isotopicdist(Info.Formula, ...
                   'Resolution', fwhh);
        end     
        %%%%%%
        mappedData = mapToMzAxis(mzVector, DF);
        % scale most abundant peak to 1
        if ~isempty(mappedData)
            mappedData = mappedData ./ max(mappedData);
            combindedPattern = combindedPattern + mappedData; 
            Q(j,:)=mappedData;
            insideInterval(j) = true;
        end        
    end
    peptidQ{sheet}=Q(insideInterval,:);
    % save pattern
    %taskName = [taskName, num2str(sheet), '.mat'];
    %save(taskName, 'combindedPattern')
    % plotting
    figure
    plot(mzVector, combindedPattern)
    disp(['Protein ' nameProteins{sheet} ' processed'])
end
peptidPattern = MSPatternCell(peptidQ,mzVector,nameProteins);