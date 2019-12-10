function [labelsAfterThr, optimalThr] = MSScore2LabelViaThreshold( scores2Transform, scores2TrainThr, groundTruth )
% This helper function transforms the MSScoreData scores2Transform to
% MSLabelData by applying thresholds that are optimal for the MSScoreData
% scores2TrainThr.
% Multi-column data are supported as well as {1xk} cell arrays of
% MSScoreData as long as k is the same for both inputs. This is useful when
% using Cross-Validation and will be used in the .performance method in
% MSClassificationResults. The output will still be {1xk}. Use .join from
% MSPredictionData to combine the cells if desired.
%
% So far only supported for the two-class scenario, i.e. a groundtruth with 
% 2 labels.

c = assertInputs( scores2Transform, scores2TrainThr, groundTruth );
saveROCCurvesFlag = false;
if ~c
  scores2Transform = packInCell( scores2Transform );
  scores2TrainThr = packInCell( scores2TrainThr );
end
numCells = numel( scores2Transform );
numPredictions = scores2Transform{1}.numPredictions;

funThrClass1Below = @(X,thr) -(X(:,1)<=thr)+2;
funThrClass1Above = @(X,thr) -(X(:,1)>=thr)+2;
optimalThr = zeros( numCells, numPredictions );

labelsAfterThr = cell(1,numCells);

for iCell = 1:numCells
  trainThrScoresSummary_iCell = MSScoreValidationSummary( scores2TrainThr{iCell}, groundTruth, saveROCCurvesFlag );
  thr_iCell = trainThrScoresSummary_iCell.optThr;
  optimalThr( iCell, : ) = thr_iCell;
  thrDirection_iCell = trainThrScoresSummary_iCell.thrDirection;
  
  scores2Transform_iCellSplitted = scores2Transform{iCell}.splitScores;
  
  for iPrediction = 1:numPredictions
    if thrDirection_iCell( iPrediction )
      labelsAfterThr_iPrediction = scores2Transform_iCellSplitted{iPrediction}.score2Label( @(X)funThrClass1Above(X,thr_iCell(iPrediction)) );
    else
      labelsAfterThr_iPrediction = scores2Transform_iCellSplitted{iPrediction}.score2Label( @(X)funThrClass1Below(X,thr_iCell(iPrediction)) );
    end
    maskNonZero = scores2Transform_iCellSplitted{iPrediction}.anyPredictionsMask;
    labelsAfterThr_iPrediction.data( ~maskNonZero, : ) = 0;
    if iPrediction == 1
      labelsAfterThr{iCell} = labelsAfterThr_iPrediction;
    else
      labelsAfterThr{iCell}.data(:,iPrediction) = labelsAfterThr_iPrediction.data;
    end
  end 
end

end

function output = packInCell( input )
  output = cell(1,1);
  output{1} = input; 
end

% Assert Input

function c = assertInputs( scores1, scores2, groundTruth )
  assertLabels( groundTruth );
  assert2Classes( groundTruth );
  c = false;
  [~,c1] = assertScoresOrScoresCA( scores1 );
  [~,c2] = assertScoresOrScoresCA( scores2 );
  if (c1 && ~c2) || (c2 && ~c1)
    error('Input type of scores not consistently scores or cell array of scores');
  end
  if c1 && c2
    c = true;
    assertMatchingNumels( scores1, scores2 );
    assertMatchingNumItems( scores1{1}, scores2{1}, groundTruth );
  else
    assertMatchingNumItems( scores1, scores2, groundTruth );
  end
end

function [e,c] = assertScoresOrScoresCA( scores )
  c = false;
  e = isa( scores, 'MSScoreData' );
  if ~e
    e = iscell( scores );
    if e
      c = true;
      for i = 1:numel(scores)
        if ~isa( scores{i}, 'MSScoreData' )
          e = false;
          break;
        end
      end
    end
  end
  handleErrors( e, 'Is not an MSScoreData or Cell Array of MSScoreData' );
end

function e = assertMatchingNumels( CAscores1, CAscores2 )
  e = numel(CAscores1) == numel(CAscores2);
  handleErrors( e, 'CAArray size does not match' );
end

function e = assertLabels( labels )
  e = isa( labels, 'MSLabelData' );
  handleErrors( e, 'Is not an MSLabelData' );
end

function e = assert2Classes( labels )
  e = labels.numLabels == 2;
  handleErrors( e, 'Groundtruth object does not have 2 labels' );
end

function e = assertMatchingNumItems( varargin )
  n = nargin;
  numItemsVec = zeros(1,n);
  for k = 1:n
    numItemsVec(k) = varargin{k}.numItems;
  end
  e = ~any( diff(numItemsVec) );
  handleErrors( e, 'NumItems do not match' );
end

function handleErrors( e, msg )
  if ~e
    error(msg);
  end
end