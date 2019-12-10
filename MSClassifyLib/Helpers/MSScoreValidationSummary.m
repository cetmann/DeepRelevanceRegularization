classdef MSScoreValidationSummary
  
  properties
    labels
    ROCAUC
    optThr = [];
    thrDirection = []; % vector of 0s and 1s, 0 == 'Class1belowThr'/ 1 == 'Class1aboveThr'
    optBalAcc = [];
    numTrainItems
    classWiseScoreMeans
    classWiseScoreQuartiles
    R = [];
  end
  
  methods
    function obj = MSScoreValidationSummary( scoreData, groundTruth, saveROCCurvesFlag )    
      obj.labels = scoreData.labels;
      R = scoreData.ROCCurve( groundTruth );
      if saveROCCurvesFlag
        obj.R = R;
      end
      obj.ROCAUC = R.computeAUC( false );
      
      if groundTruth.numLabels == 2
      
        obj.optThr = R.computeOptThr;

        % Direction of threshold - this information is related to the
        % positive class to be defined in MSComputeROCCurve, and there 
        % should exist a better way to determine it or fix it ?!
        mask1 = scoreData.anyPredictionsMask & groundTruth.data == 1;
        mask2 = scoreData.anyPredictionsMask & groundTruth.data == 2;
        obj.thrDirection = median( scoreData.data(mask1,1:2:end), 1 ) > median( scoreData.data(mask2,1:2:end), 1 );
        obj.optBalAcc = R.computeOptBalAcc;
      
      end
      
      obj.numTrainItems = sum(scoreData.anyPredictionsMask);
      idxLabel1 = groundTruth.data == 1;
      idxLabel2 = groundTruth.data == 2;
      scoresLabel1 = scoreData.data(idxLabel1,1:2:end);
      scoresLabel2 = scoreData.data(idxLabel2,1:2:end);
      cWSM1 = mean( scoresLabel1, 1 ); 
      cWSM2 = mean( scoresLabel2, 1 );
      p = [0.25 0.5 0.75];
      cWSQ1 = prctile( scoresLabel1, p, 1 );
      cWSQ2 = prctile( scoresLabel2, p, 1 );
      obj.classWiseScoreMeans = [cWSM1;cWSM2];
      obj.classWiseScoreQuartiles = [cWSQ1; cWSQ2];
    end
  end
  
end
