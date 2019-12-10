function dataPartition=MSGenerateDataPartition(maldiData, classesNames,...
     classesPattern, sampleSetsPattern, samplesPattern, individualsPattern,...
                       cvEntities, balanced)
% function used to generate a data partition from input patterns. The
% classes label, representing the ground truth is always generated.
% Other labels like samples (cores), sampleSets (TMAs), and individuals
% (patients) are built depending on whether the corresponding patterns are
% given. The label cvEntities is a reference to one of the generated labels
% or empty (user's choice) and is used for generating the cv-segments
% for the crossvalidation task.
%
%INPUT
%   maldiData: MSMaldiData object containing the training data
%   classesNames: cell array of strings with dimension number of classes
%                 containing the labels of the classes.
%                 e.g. classesNames = {'PancreasTumor','LungTumor'}
%   classesPattern: cell array of strings with dimension number of classes
%                 containing the regular expressions for class extraction
%                 e.g. classesPattern = {'AD_','T_'}  
%   sampleSetsPattern: cell array of strings containing the regular 
%                 expressions for sampleSets (TMAs) extraction
%                 e.g. sampleSetsPattern={'/L8_1_Mark_Kriegsmann$'; ...
%                                         '/L1_1_MarkKriegsmann$'; ...
%                                         '/D7_panc_3A$';...
%                                         '/TMA_D04_panc$'};
%   samplesPattern: String containing the regular 
%                 expressions for samples (cores) extraction
%                 e.g. regionPattern =...
%                    '(?:Kriegsmann/\d+|L\d_\d/\d+|_panc/\d+|_panc_3A/\d+)'
%   individualsPattern: cell array of strings containing the regular 
%                 expressions for individuals (patients) extraction. This
%                 regular expression should match the last part of the
%                 annotation name, right after reading the last character
%                 that is part of the name. Here the 'look-behind' pattern
%                 '(?<=test)expr' could be useful.
%                 e.g. individualsPattern=...
%                                       {'(?<=AD_\d+)_.*','(?<=T_\d+)\D.*'}
%   cvEntities:   Either empty or one of these three strings:"samples", 
%                 "sampleSets", "individuals". Indicates which label is
%                 used to generate the cv-segments for crossvalidation
%   balanced:     logical scalar which indicates whether the different
%                 classes should be balanced (same number of samples)
% OUTPUT
%   dataPartition: MSDataPartition object

%-------------------- input validation ------------------------------------
narginchk(3,8)
if nargin<8
    balanced=true;
elseif ~islogical(balanced)||~isscalar(balanced)
        error('input argument "balanced" must be a logical scalar')
end
if nargin < 7
    cvEntities=[];
end
if nargin < 6
    individualsPattern=[];
else
    validate(individualsPattern);
end
if nargin < 5
    samplesPattern=[];
else
    if~ischar(samplesPattern)
        error('samplesPattern must be a string')
    end
end
if nargin < 4
    sampleSetsPattern=[];
else
    validate(sampleSetsPattern);
end
if ~isa(maldiData,'MSMaldiData')
    error('maldiData must be an instance of MSMaldiData')
else
    maldiData.assert();
end
validate(classesNames);
validate(classesPattern);
if length(classesNames)~=length(classesPattern)
    error(['The length of classesNames and classesPattern cellarrays'...
          'must be the same'])
end

%----------------------- initialize data partition ------------------------
dataPartition=MSDataPartition();
%initialize classes pattern
classPatternString=classesPattern{1};
for i=2:length(classesPattern)
    classPatternString=strcat(classPatternString,'|',classesPattern{i});
end
%----------------------- data partition classes ---------------------------
dataPartition.classes = maldiData.partition(classesNames, classesPattern); 
% if balanced data option is active
if balanced
    %get labels
    usedLabels=dataPartition.classes.usedLabels(dataPartition.classes.usedLabels>0);
    %build a hash for counting labels
    hash=zeros(1,usedLabels(end));
    %count number of occurences per class
    classes=dataPartition.classes.data(dataPartition.classes.data>0);
    for i=1:length(classes)
        hash(classes(i))=hash(classes(i))+1;
    end
    %find least number of members per class
    minim=min(hash(usedLabels));  
    %even classes
    for i=1:length(usedLabels)
        difference=hash(usedLabels(i))-minim;
        if difference>0
            indices = find(dataPartition.classes.data == usedLabels(i));
            randIndices = randperm(size(indices,1),difference);
            dataPartition.classes.data(indices(randIndices)) = 0;
        end
    end
end
%------ data partition annotated samples (stored in other struct) ---------
[tumorIndices, tumorNames] = maldiData.annotations.find(classPatternString);
% Create label for annotated samples (cores)
dataPartition.other.annotatedSamples= maldiData.partition(tumorNames, ...
                                                   num2cell(tumorIndices));
%------------------------ data partition sampleSets -----------------------
if ~isempty(sampleSetsPattern)
    % Extract TMA annotations
    sampleSetsIndices = [];
    sampleSetsNames = [];
    for i=1:size(sampleSetsPattern,1)
        [sampleSetsIndicesNew, sampleSetsNamesNew] = maldiData.annotations.find(...
            sampleSetsPattern{i});
        sampleSetsIndices = [sampleSetsIndices , sampleSetsIndicesNew];
        sampleSetsNames = [sampleSetsNames, sampleSetsNamesNew];
    end
    % Create label for samples (cores)
    dataPartition.sampleSets = maldiData.partition(sampleSetsNames,...
        num2cell(sampleSetsIndices)); 
end
%------------------------ data partition samples --------------------------
if ~isempty(samplesPattern)
    [regionIndices, regionNames] = maldiData.annotations.find(samplesPattern);
    % Create label for samples (cores)
    dataPartition.samples = maldiData.partition(regionNames,num2cell(regionIndices));
end
%------------------------ data partition individuals ----------------------
if ~isempty(individualsPattern)
    % Each string in individualsPattern is a regular expression returning
    % string behind the patient name
    if ~isempty(individualsPattern)
        patientNames = regexprep(tumorNames,individualsPattern{1},'');
        for i=2:length(individualsPattern)
            patientNames = regexprep(patientNames,individualsPattern{i},'');%Lung TMAs
        end
    end

    % Find unique patient names and the map from core label to patient label
    [patientNames, ~, patientLabelMap] = unique(patientNames);
    % Create label vector for patient labels
    patientLabels = zeros(size(dataPartition.other.annotatedSamples.data));
    nonZeroLabels = (dataPartition.other.annotatedSamples.data > 0);
    patientLabels(nonZeroLabels) = patientLabelMap(...
    dataPartition.other.annotatedSamples.data(nonZeroLabels));
    dataPartition.individuals = MSLabelData(patientNames,...
        patientLabels, dataPartition.other.annotatedSamples);
    % Map patient ids to entire cores (trainingPartition.samples) instead
    % of confined id to annotated tumor region
    aggregatedIndividuals = dataPartition.individuals.aggregate(...
       dataPartition.samples);
    dataPartition.individuals = aggregatedIndividuals.aggregateInverse( ...
       dataPartition.samples);
end
%-------------------------- data partition cvEntities ---------------------
if ~isempty(cvEntities)
    try
        dataPartition.cvEntities=dataPartition.(cvEntities);
    catch
        error(['cvEntities must be one of these three strings:'...
              ' "samples", "sampleSets", "individuals"'])
    end
end
end

function validate(pattern)
if ~(isempty(pattern)||(iscellstr(pattern)&& isvector(pattern)))
    error('The pattern should be a cell array of strings')
end
end

