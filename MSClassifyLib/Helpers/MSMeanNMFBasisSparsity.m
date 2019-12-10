function meanSparcity=MSMeanNMFBasisSparsity(nmfMap, numBasis, plotFlag)
narginchk(1,3)
nBasis=size(nmfMap.basis,1);
plotFlag=(nargin<3)||(islogical(plotFlag)&&isscalar(plotFlag)&&plotFlag);
if nargin<2
    numBasis=nBasis;
else
    numBasis=min(nBasis,numBasis);
end
nmfMap.sortBasisVectors('l1l2');
meanSparcity=mean(nmfMap.featureMapInfo.l1l2Norms(1:numBasis));
if plotFlag
    plot(nmfMap.featureMapInfo.l1l2Norms),hold on
end