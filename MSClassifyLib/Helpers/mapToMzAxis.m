function [mappedData] = mapToMzAxis(mzVector, DF)
% Maps isotopic pattern to given mz-vector
% input: mzVector target mz-axis
%        DF density function of isopic pattern

    minMz = min(mzVector) - 1;
    maxMz = max(mzVector) + 1;
    minIsoMz = min(DF(:,1));
    maxIsoMz = max(DF(:,1));
    
    % Return empty mappedData if mzVector does not contain pattern interval
    % completely
    if minMz>minIsoMz||maxMz<maxIsoMz
        mappedData=[];
        return
    end
    
    % create data left/right of pattern with zeros with a resolution given
    % by the targeted mz-axis
    resMz = diff(mzVector);  % take resolution from target mzVector
    isoMzVectorL = minMz: resMz(1,1) : minIsoMz - resMz(1,1) ;
    isoMzVectorR = maxIsoMz + resMz(1,1): resMz(1,1) : maxMz;
    
    % concatenated zero-regions on the left/right with isotopic pattern
    isoDataL = zeros(2, size(isoMzVectorL,2));
    isoDataL(1,:) = isoMzVectorL;
    isoDataR = zeros(2, size(isoMzVectorR,2));
    isoDataR(1,:) = isoMzVectorR;
    isoData = [isoDataL , DF' , isoDataR];
    mappedData = MSResampleData(isoData(2,:),isoData(1,:),mzVector);
end