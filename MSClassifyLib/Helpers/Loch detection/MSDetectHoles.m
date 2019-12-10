function mask=MSDetectHoles(maldiData, holePeakCriteria, tissuePeakCriteria, logicFlag, criteria, postProcFlag)

% Return a mask (logical array) which represents the classification of holes (true) and
% tissue (false)
% mask=MSDetectHoles(maldiData, holePeakCriteria, tissuePeakCriteria, logicFlag, criteria)
%     INPUT  
%     maldiData: Object of type MSMaldiData
%     holePeakCriteria (tissuePeakCriteria): struct with fields:
%       -peaks: array with the relevant peaks only mostly present in hole
%       (tissue) spectra.
%       -radius: positive scalar specifying the size of the neighbourhood
%       around peaks to be considered
%       -threshold: positive scalar specifying the minimum (maximum) value
%       of the sum or norm for including the spectra in the mask
%     logicFlag: string '&&' or '||' to decide whether the union or
%     intersection of the spectra found with each PeakCriteria is performed
%     criteria: string 'norm', 'sum', 'norm_rel_mean' or 'sum_rel_mean' 
%     determining which criteria will be used to choose the spectra    
%
%     OUTPUT
%     mask:    Logical array representing hole (true) and tissue (false) cells
%
%     Selection of peaks among those shown bellow is recommended for hole
%     detection:
%     holePeaks=[842.5 1012]
%     tissuePeaks=[868.5 852.5 840.5 836.5 898.5 1198.5 1138.5 1105.8 976.5]

% input validation
narginchk(5,6);
if nargin < 6
    postProcFlag=1;    
elseif ~(isscalar(postProcFlag) && isnumeric(postProcFlag))
    error('The post processing flag must be a numeric scalar quantity')
end
if nargin < 5
  criteria = 'norm';
else
criteriaList={'norm_rel_mean', 'sum_rel_mean', 'norm', 'sum'};
if ~any(strcmp(criteriaList, criteria))
    error('The selection criteria is not valid')
end
end

if nargin < 4
   logicFlag='&&';
elseif ~(strcmp(logicFlag,'&&')||strcmp(logicFlag,'||'))
    error('LogicFlag must be one of the two strings: "&&" or "||"');
end

maldiData.assert
verify(holePeakCriteria);
verify(tissuePeakCriteria);

% mzValues in intervals centered in each mzValue in peaks 
holePeaksSpect=mzValues(maldiData.mzVector,holePeakCriteria.peaks,holePeakCriteria.radius);
tissuePeaksSpect=mzValues(maldiData.mzVector, tissuePeakCriteria.peaks, tissuePeakCriteria.radius);

auxcount=~isempty(holePeaksSpect)+2*~isempty(tissuePeaksSpect);%0->both empty, 1->tissue empty, 2->hole empty, 3->non empty

% decide whether the threshold is applied to the original spectra or
% relative to the mean
if strcmp(criteria,'norm_rel_mean')||strcmp(criteria,'sum_rel_mean')
    f=@findSpectrRelative2Mean;
else
    f=@findSpectr;
end    
%obtain the mask based on the threshold criteria
switch auxcount
    case 0
        mask=false(size(maldiData.data,1));
        return
    case 1
        indexes=f(maldiData.data(:,holePeaksSpect), holePeakCriteria.threshold, criteria, '>');
    case 2
        indexes=f(maldiData.data(:,tissuePeaksSpect), tissuePeakCriteria.threshold, criteria, '<');
    case 3
        indexesH=f(maldiData.data(:,holePeaksSpect), holePeakCriteria.threshold, criteria, '>');
        indexesT=f(maldiData.data(:,tissuePeaksSpect), tissuePeakCriteria.threshold, criteria, '<');
        switch logicFlag
            case '&&'
                indexes=indexesH & indexesT;
            case '||'
                indexes=indexesH | indexesT;
        end
end
% Eliminate isolated hole-spectra

%computing closure of the image representing the two classes (hole, tissue) 
if postProcFlag    
    I=maldiData.positions.encube(indexes'+0);
    I = medfilt2(I);
    %seMatrix=[0 1 0;1 1 1; 0 1 0];
    %seMatrix=[1 1;1 1];
    seMatrix=1;
    se = strel('arbitrary',seMatrix);
    %se = strel('disk',1,0) 
    closeI = imclose(I,se);
    %applying median filtering to eliminate 'salt & pepper' effect
    %K = medfilt2(closeI);
    K=closeI;
    %eliminate conex components with few pixels
    K=bwareaopen(K,15);
    %the mask is obtained, where 1 corresponds to holes and 0 to tissue
    mask=logical(K(maldiData.positions.reverseIndex));
else
    mask=indexes';
end
end

function y=mzValues(mzVector, mzPeaks, radius)
N=size(mzPeaks,2);
for i=1:N
    % find the interval of rows in data corresponding to the m/z value interval
    centerMZ=mzPeaks(i);
    iLeft=MSFindMZIndex(mzVector, centerMZ-radius);
    iRight=MSFindMZIndex(mzVector, centerMZ+radius);
    y(i)={iLeft:iRight};
end
if exist('y','var')
    y=unique(cell2mat(y));
else
    y=[];%in case the vector is empty (meaning that the mzPeak is empty)
end
end

function verify(peakCriteria)
% Checks whether the peakCriteria parameter is an appropriate struct and
% sends an error message in case one of the conditions is violated
%   verify(peakCriteria)
%   INPUT
%   peakCriteria: : struct with fields
%       -peaks: array with positive entries
%       -radius: positive scalar 
%       -threshold: positive scalar 
if ~(isstruct(peakCriteria) && isfield(peakCriteria,'peaks')...
        && isfield(peakCriteria, 'threshold') && isfield(peakCriteria, 'radius')...
        && isnumeric(peakCriteria.peaks) && isnumeric(peakCriteria.threshold)...
        && isnumeric(peakCriteria.radius) && isscalar(peakCriteria.radius) && isscalar(peakCriteria.threshold)...
        && peakCriteria.radius>=0 && peakCriteria.threshold>=0) 
    error('Expected a structure with numeric positive fields \"peaks\", \"threshold\" and \"radius\" ')
end
end

function indexes=findSpectrRelative2Mean(data, threshold, criteria, compareFlag)
% Return a mask for the the spectra with high (or low) intensity (when compared to
% the mean) according to the specified threshold and criteria 
%   indexes=findSpectrRelative2Mean(data, threshold, criteria, compareFlag)
%     INPUT
%     dataHole: Two-dimensional data (Spectra vs m/z) related to hole
%     dataTissue: Two-dimensional data (Spectra vs m/z) related to tissue
%     threshold: scalar; the spectra above threshold*mean_spectra is chosen
%     criteria: String 'norm' or 'sum' determining which criteria will be
%     used to choose the spectra
%     compareFlag: string '>' or '<' indicating whether spectra with small
%     or large norm or sum are wanted
%
%     OUTPUT
%     indexes:    Indexes of the spectra with relevant intensity

if ~(ismatrix(data)&&isnumeric(data)&&~isempty(data))
    error('Data matrix must be a non-empty, numeric 2D-array')
end

% processing dataHole
[nrows,~]=size(data);

meanData=mean(data,1);% mean spectrum for the interval of interest
targetUnbiasSpect=data-repmat(meanData,nrows,1); %substracting the mean from spectrum
if compareFlag=='>'
    %numel=sum(targetUnbiasSpect>0);
    targetUnbiasSpect(targetUnbiasSpect<0)=0; %eliminating negative entries (only peaks above mean are important)
else
    %numel=sum(targetUnbiasSpect<0);
    targetUnbiasSpect(targetUnbiasSpect>0)=0; %eliminating positive entries (only peaks bellow mean are important)
    targetUnbiasSpect=-targetUnbiasSpect;
end

switch criteria
    case 'norm_rel_mean'
        spect=(sqrt(sum((targetUnbiasSpect).^2,2)))'; %computing the norm for each row 
    case 'sum_rel_mean'
        spect=sum(targetUnbiasSpect,2)'; %computing the elements' sum for each row
end
% threshold means the percentage of highest(lowest) intensity
threshold=min(threshold,1);
[~, I]=sort(spect,'descend');
indexes=false(size(spect));
indexes(I(1:floor(threshold*end)))=true;

% threshold means an intensity fraction of the highest intensity spectra
%indexes=(spect/max(spect)>threshold); 
end

function indexes=findSpectr(data, threshold, criteria, compareFlag)
% Return a mask for the the spectra with high (or low) intensity (when compared to
% the mean) according to the specified threshold and criteria 
%   indexes=FindSpectrRelative2Mean(data, threshold, criteria, compareFlag)
%     INPUT
%     dataHole: Two-dimensional data (Spectra vs m/z) related to hole
%     dataTissue: Two-dimensional data (Spectra vs m/z) related to tissue
%     threshold: scalar; the spectra above threshold*mean_spectra is chosen
%     criteria: String 'norm' or 'sum' determining which criteria will be
%     used to choose the spectra     
%     compareFlag: string '>' or '<' indicating whether spectra with small
%     or large norm or sum are wanted
%
%     OUTPUT
%     indexes:    Indexes of the spectra with high intensity

if ~(ismatrix(data)&&isnumeric(data)&&~isempty(data))
    error('Data matrix must be a non-empty, numeric 2D-array')
end
switch criteria
    case 'norm'
        spect=(sqrt(sum((data).^2,2)))'; %computing the norm for each row        
    case 'sum'
        spect=sum(data,2)'; %computing the elements' sum for each row
end


% threshold means the percentage of highest(lowest) spectra intensity
threshold=min(threshold,1);
switch compareFlag %selecting significant indexes
    case '>'
        [~, I]=sort(spect,'descend');
    case '<'
        [~, I]=sort(spect,'ascend');
end      

indexes=false(size(spect));
indexes(I(1:floor(threshold*end)))=true;

% threshold means an intensity fraction of the highest intensity
% spect=spect/max(spect);
% switch compareFlag %selecting significant indexes
%     case '>'
%         indexes=(spect>threshold); 
%     case '<'
%         indexes=(spect<threshold);
% end       
end