function [NLL,l1,l2,l1l2,logitSens,logitSqSens,l1l2LogitSens,probSens,lossSens] = MSReadNNLogFile(nnModel)
    nfName=nnModel.neuralNetFile;
    [pathstr,name,~] = fileparts(nfName);
    logFile = [pathstr name '_log.mat'];
    NLL = h5read(logFile,'/NLL');
    l1 = h5read(logFile,'/l1');
    l2 = h5read(logFile,'/l2');
    l1l2 = l1./(sqrt(l2+10*eps));
    logitSens = h5read(logFile,'/logitSensitivity');
    logitSqSens = h5read(logFile,'/logitSqSensitivity');
    l1l2LogitSens = logitSens./(sqrt(logitSqSens +10*eps));
    probSens = h5read(logFile,'/probabilitySensitivity');
    lossSens = h5read(logFile,'/NLLSensitivity');
end
