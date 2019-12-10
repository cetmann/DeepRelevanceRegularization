%% Set file paths %%
addpath(genpath('/path/to/MSClassifyLib'))
load('/path/to/AP43_B12_O123_HB_TR_afx_newAnnotation.mat')

%% Prepare data
data.setNormalization('tic')
data.reduceMzRange([700,2000])

P = partition;
D = data;


nFolds=5;
trainCenter = {'TR','HB'};
testCenter = {'HB','TR'};


%%
onlyFineAnnotations = (P.sampleSets.data<9);
onlyWithPatientNumbers = (P.individuals.data~=0);

P.sampleSets.data = P.sampleSets.data .* (P.individuals.data~=0);
CV=P.sampleSets.kFold(nFolds,P.individuals);
labs = MSLabelData(...
    {'HB','TR'},...
    mod(P.sampleSets.data,2) + 2.*(mod(P.sampleSets.data,2)==0) .* (P.sampleSets.data~=0),...
    D); % odd numbers correspond to HB, even correspond to TR




parameters = [num2cell(5) num2cell(10:10:100) num2cell(125:25:200)];
nParameters = length(parameters);
nCenters = length(trainCenter);

ind = 1:nFolds;

experiment = struct();
experiment.parameters = parameters;
experiment.center = cell(2,1);
experiment.predictionLabels = zeros(size(P.classes.data));
experiment.classes = MSLabelData(...
    {'B','O'},...
    P.classes.data .* onlyFineAnnotations .* onlyWithPatientNumbers,...
    D);



% We use the "i=1;while i<N;i++" scheme over for loops to be able to
% restart the CV scheme at any point.

%% Run calculations %%

iCenter = 1;
while iCenter<=2
    experiment.center{iCenter} = struct();
    experiment.center{iCenter}.trainCenter = trainCenter{iCenter};
    experiment.center{iCenter}.testCenter = testCenter{iCenter};
    experiment.center{iCenter}.testFold = cell(nFolds,1);
    experiment.center{iCenter}.predictionLabels = zeros(size(experiment.classes.data));
    
    iTestFold = 1;
    while iTestFold <= nFolds
        experiment.center{iCenter}.testFold{iTestFold} = struct();
        testLabels = experiment.classes.data .* (labs.data==iCenter) .* (CV.data==iTestFold) .* onlyFineAnnotations .* onlyWithPatientNumbers;
        experiment.center{iCenter}.testFold{iTestFold}.testLabels = MSLabelData({'B','O'},testLabels,D);
        experiment.center{iCenter}.testFold{iTestFold}.P = cell(nParameters,1);
        experiment.center{iCenter}.testFold{iTestFold}.valIndices = ind(1:end ~=iTestFold);
        experiment.center{iCenter}.testFold{iTestFold}.validationScore = 0;
        
        iParameter = 1;
        while iParameter <= nParameters
            experiment.center{iCenter}.testFold{iTestFold}.P{iParameter} = struct();
            experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.valFold = cell(nFolds-1,1);
            experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.predictionLabels = zeros(size(experiment.classes.data));
            
            iValFold = 1;
            while iValFold < (nFolds-1)
                
                valSampleSetIndex = ...
                    experiment.center{iCenter}.testFold{iTestFold}.valIndices(iValFold);
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.valFold{iValFold}.valSampleSetIndex = valSampleSetIndex;
                
                valLabels = experiment.classes.data .* ...
                    (labs.data==iCenter) .* ...
                    (CV.data==valSampleSetIndex) .* ...
                    onlyFineAnnotations .* ...
                    onlyWithPatientNumbers;
                
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.valFold{iValFold}.valLabels = ...
                    MSLabelData({'B','O'},valLabels,D);
                
                
                trainLabels = experiment.classes.data .* ...
                    (labs.data~=iCenter) .* ...
                    (CV.data~=valSampleSetIndex) .* ...
                    (CV.data~=iTestFold) .* ...
                    onlyFineAnnotations .* ...
                    onlyWithPatientNumbers;
                
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.valFold{iValFold}.trainLabels = ...
                    MSLabelData({'B','O'},trainLabels,D);
                
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.valFold{iValFold}.trainSampleSetIndices = ...
                    ind(1:end~=valSampleSetIndex & 1:end~=iTestFold);
                
                parameter = parameters{iParameter};
                
                
                roc = MSROCPicking(parameter);
                rocExtractor = roc.createMap(D, trainLabels);
                
                rocTrainFeatures = rocExtractor.map(D, trainLabels);
                rocValFeatures = rocExtractor.map(D, valLabels~=0);
                
                
                
                classifier = MSLDAClassifier();
                
                
                
                model = classifier.trainModel(rocTrainFeatures, trainLabels);
                predictionLabels = model.classify(rocValFeatures, valLabels~=0);
                
                
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.valFold{iValFold}.predictionLabels = ...
                    predictionLabels;
                
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.predictionLabels = ...
                    experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.predictionLabels + 1.0 * predictionLabels.data;
                
                iValFold = iValFold + 1;
                
                save(['ROCLDA_BO_center_',num2str(iCenter),'.mat'],'experiment','iCenter','iParameter','iValFold','iTestFold','-v7.3')
            end
            experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.predictionLabels = ...
                MSLabelData({'B','O'},...
                experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.predictionLabels,...
                D);
            
            score = evaluate(experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.predictionLabels,P);
            experiment.center{iCenter}.testFold{iTestFold}.P{iParameter}.score = score;
            
            if iParameter == 1
                experiment.center{iCenter}.testFold{iTestFold}.bestParameter = parameters{iParameter};
                experiment.center{iCenter}.testFold{iTestFold}.bestParameterIndex = iParameter;
            else
                experiment.center{iCenter}.testFold{iTestFold}.bestParameter = parameters{iParameter-1};
                experiment.center{iCenter}.testFold{iTestFold}.bestParameterIndex = iParameter-1;
                
            end
            iParameter = iParameter + 1;
        end
        
        
        
        
        
        trainLabels = experiment.classes.data .* ...
            (labs.data~=iCenter) .* ...
            (CV.data~=iTestFold) .* ...
            onlyFineAnnotations .* ...
            onlyWithPatientNumbers;
        experiment.center{iCenter}.testFold{iTestFold}.trainLabels = ...
            MSLabelData({'B','O'},...
            trainLabels,...
            D);
        
        roc = MSROCPicking(parameter);
        rocExtractor = roc.createMap(D, trainLabels);
        
        rocTrainFeatures = rocExtractor.map(D, trainLabels);
        rocTestFeatures = rocExtractor.map(D, experiment.center{iCenter}.testFold{iTestFold}.testLabels);
        
        
        
        classifier = MSLDAClassifier();
        
        
        
        model = classifier.trainModel(rocTrainFeatures, trainLabels);
        
        
        
        
        predictionLabels = model.classify(rocTestFeatures, experiment.center{iCenter}.testFold{iTestFold}.testLabels);
        
        
        
        
        experiment.center{iCenter}.testFold{iTestFold}.predictionLabels = ...
            predictionLabels;
        
        experiment.center{iCenter}.predictionLabels = ...
            experiment.center{iCenter}.predictionLabels + predictionLabels.data;
        
        
        
        
        
        
        iTestFold = iTestFold + 1;
        save(['ROCLDA_BO_center_',num2str(iCenter),'.mat'],'experiment','iCenter','iParameter','iValFold','iTestFold','-v7.3')
    end
    experiment.center{iCenter}.predictionLabels = ...
        MSLabelData({'B','O'},...
        experiment.center{iCenter}.predictionLabels,...
        D);
    experiment.center{iCenter}.score = evaluate(...
        experiment.center{iCenter}.predictionLabels,P);
    
    
    experiment.predictionLabels = ...
        experiment.predictionLabels + experiment.center{iCenter}.predictionLabels.data;
    
    
    iCenter = iCenter + 1;
end
experiment.predictionLabels = ...
    MSLabelData({'B','O'},...
    experiment.predictionLabels,...
    D);
experiment.score = evaluate(experiment.predictionLabels,P);


function score = evaluate(results,P)
evaluationMetric = 'avgSensitivity';
[~,confusionMatrix] = results.validate(...
    P.classes.data .* (results.data~=0));
confusionMeasure = MSConfusionMeasure(confusionMatrix);
score = confusionMeasure.(evaluationMetric);
end
