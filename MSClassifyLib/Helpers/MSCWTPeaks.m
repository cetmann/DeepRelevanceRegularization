function [peakVec,sNRVec,ridgeLengthVec] =  MSCWTPeaks(spectrum,varargin)
% This function calculates the indices of the peaks in a given data set
% using continuous wavelet transformation.
%
% Parameters:
%
% MANDATORY:
% spectrum - A series of measurements. The spectrum can be of type single
% or double
%
% OPTIONAL:
% Scales : 
%
% The scale values for the continuous wavelet transformation
%
% Wavelet: 
%
% The wavelet used in the CWT
%
% FilterLenght: 
%
% The length of the filter which is used for smoothing the
% CWT coefficients matrix
%
% Method:  
% 
% The peak position selection method. Either the peak position is
% determined  by finding the index of the maximal CWT coefficient
% ('maxCoef') or the peak index correspondes to the initial index of the
% peaks ridge in the CWT matrix ('initialIndex')
% 
% NoiseThreshold:
%
% If a peaks noise is lower than this threshold, a default values is used.
%
% NoiseWindowSize:
%
% The size of the sliding window, when calculating the peaks noise.
%
% SNRThreshold:
%
% Peaks with Signal-to-noise-ratio lower than this threshold are sorted
% out.
%
% LengthThreshold:
%
% Peaks, whose length is lower than this threshold, are sorted out.
%
% MinAmpScale and MaxAmpScale
%
% When sorting out unwanted peaks, these values determine the minimal /
% maximal scale in which maximal CWT coefficients are allowed. If they do
% not lie in this interval, these peaks are sorted out.
%
% Plot:
%
% If 'Plot' is set to 'true', a plot of the CWT coefficients is given.
%
% OUTPUT:
% peakVec - A vector containing the indices of the peaks in the input spectrum
%
% The optional parameters can either be passed by string reference, e.g.
% cWTPeaks(spectrum,'Scales', 1:20,'Method','initialIndex',…) or in the order, in which they appear in
% the source code.
%
%Reference: Du P., Kibbe W.A. and Lin S.M. (2006), Improved Peak Detection in Mass Spectrum by Incorporating
%Continuous Wavelet Transform-based Pattern Matching, Bioinformatics,
%Volume 22, Issue 17. See also:
%http://basic.northwestern.edu/publications/peakdetection/

%% Parse Input
%Set default values
DEFAULT_SCALES = 1:13;
DEFAULT_NOISE_THRESHOLD = mean(spectrum);
DEFAULT_SNR_THRESHOLD = 1;
DEFAULT_FILTER_LENGTH = 1;
DEFAULT_WAVELET = 'mexh';
DEFAULT_LENGTH_THRESHOLD = 1;
DEFAULT_MIN_AMP_SCALE = 1;
DEFAULT_MAX_AMP_SCALE = numel(DEFAULT_SCALES);
DEFAULT_METHOD = 'maxCoef';
DEFAULT_NOISE_WINDOW_SIZE = 50;

% Parse input
p = inputParser;
addRequired(p,'spectrum',@isnumeric);
addOptional(p,'Scales', DEFAULT_SCALES, @isnumeric);
addOptional(p,'NoiseThreshold', DEFAULT_NOISE_THRESHOLD, @isnumeric);
addOptional(p,'SNRThreshold', DEFAULT_SNR_THRESHOLD, @isnumeric);
addOptional(p,'FilterLength', DEFAULT_FILTER_LENGTH, @isnumeric);
addOptional(p,'Wavelet', DEFAULT_WAVELET, @ischar);
addOptional(p,'LengthThreshold', DEFAULT_LENGTH_THRESHOLD, @isnumeric);
addOptional(p,'MinAmpScale',DEFAULT_MIN_AMP_SCALE,@isnumeric);
addOptional(p,'MaxAmpScale',DEFAULT_MAX_AMP_SCALE,@isnumeric);
addOptional(p,'Plot',false,@islogical);
addOptional(p,'Method',DEFAULT_METHOD,@ischar);
addOptional(p,'NoiseWindowSize',DEFAULT_NOISE_WINDOW_SIZE,@isnumeric);

parse(p,spectrum,varargin{:});

%Set parsed values
spectrum        = p.Results.spectrum;
SCALES          = p.Results.Scales;
NOISE_THRESHOLD = p.Results.NoiseThreshold;
NOISE_WINDOW_SIZE = p.Results.NoiseWindowSize;
SNR_THRESHOLD   = p.Results.SNRThreshold;
FILTER_LENGTH   = p.Results.FilterLength;
WAVELET         = p.Results.Wavelet;
LENGTH_THRESHOLD= p.Results.LengthThreshold;
MIN_AMP_SCALE   = p.Results.MinAmpScale;
if p.Results.MaxAmpScale > max(SCALES)
    MAX_AMP_SCALE  = max(SCALES);
else
    MAX_AMP_SCALE  = p.Results.MaxAmpScale;
end
METHOD          = p.Results.Method;

%% Check for correctness of input arguments
if numel(SCALES) < 2 || prod(size(SCALES) > 1) ~= 0
    error('The value of "Scales" is invalid. Please keep in mind that "Scales"  must be a vector with at least two elements.')
end

if NOISE_THRESHOLD <= 0 || numel(NOISE_THRESHOLD) > 1
       error('The value of "NoiseThreshold" is invalid. Please keep in mind that "NoiseThreshold" must be a single value bigger than zero.');
end

if SNR_THRESHOLD <= 0 || numel(SNR_THRESHOLD) > 1
       error('The value of "SNRThreshold" is invalid. Please keep in mind that "SNRThreshold" must be a single value bigger than zero.');
end

if FILTER_LENGTH <= 0 || numel(FILTER_LENGTH) > 1 || mod(FILTER_LENGTH,1) ~= 0 || mod(FILTER_LENGTH,2) == 0
       error('The value of "FilterLength" is invalid. Please keep in mind that "FilterLength" must be a single odd integer value bigger than zero.');
end

if LENGTH_THRESHOLD <= 0 || numel(LENGTH_THRESHOLD) > 1 || mod(LENGTH_THRESHOLD,1) ~= 0
       error('The value of "LengthThreshold" is invalid. Please keep in mind that "LengthThreshold" must be a single integer value bigger than zero.');
end

if MIN_AMP_SCALE < 1 || MIN_AMP_SCALE > max(SCALES) - 1 || numel(MIN_AMP_SCALE) > 1 || mod(MIN_AMP_SCALE,1) ~= 0
    error('The value of "MinAmpScale" is invalid. Please keep in mind that "MinAmpScale" must be a single integer value bigger than one and less than numel(Scales) - 1.');
end

if MAX_AMP_SCALE < 2 || numel(MAX_AMP_SCALE) > 1 || mod(MAX_AMP_SCALE,1) ~= 0
    error('The value of "MaxAmpScale" is invalid. Please keep in mind that "MaxAmpScale" must be a single integer value bigger than two and less than numel(Scales).');
end

if NOISE_WINDOW_SIZE <= 2 || numel(NOISE_WINDOW_SIZE) > 1 || mod(NOISE_WINDOW_SIZE,1) ~= 0 || NOISE_WINDOW_SIZE > numel(spectrum)/2 - 1 
       error('The value of "NoiseWindowSize" is invalid. Please keep in mind that "NoiseWindowSize" must be a single integer value bigger than two. Also values bigger than numel(spectrum)/2 - 1 are not allowed');
end

if ~(strcmpi(METHOD,'maxCoef') || strcmpi(METHOD,'initialIndex'))
   error('The value of "Method" is invalid. It must either be "maxCoef" or "initialIndex".'); 
end


%% Set-up

% Compute continuous wavelet transform coefficients and smoothen them if
% necessary
if p.Results.Plot
    coefs = cwt(spectrum,SCALES,WAVELET,'plot');
else
    coefs = cwt(spectrum,SCALES,WAVELET);
end

if FILTER_LENGTH > 1
    coefs = conv2(coefs,fspecial('average',[1 FILTER_LENGTH]),'same');
end

%Find local maxima in coefficient matrix

[maxMat,~] = localmax(coefs,[],false);

beginMaxMat = maxMat;

%% Ridge detection

%(1) Initialize ridges struct with the fields:
% nextIndex     - The peaks next index in the localmax matrix
% initialIndex  - Initial index on the first scale
% length        - The peaks length
% maxCoef       - Maximal coefficient in the CWT matrix
% maxCoefIndex  - MaxCoef's index
% maxCoefScale  - MaxCoef's scale
% SNR           - Signal to noise ratio

initialZeros = repmat({0},numel(spectrum),1);
global ridges
ridges = struct('nextIndex',[], 'initialIndex',[],...
    'SNR',[], 'length',[],... 
    'maxCoef', initialZeros, 'maxCoefIndex', [],'maxCoefScale',[]);
%Set size of ridges struct
ridges(numel(spectrum));



%Filter out duplicated ridge beginnings: If one ridge has several beginning,
%close to each other, only the beginning closest to the ridge is kept.
ridgeBeginningIndices = zeros(1,numel(spectrum));
index = 0;
for k = unique(nonzeros(maxMat(1,:))')
    appearances = find(maxMat(1,:)==k);
    [~,minDistIndex] = min(abs(k-appearances));
    index = index + 1;
    ridgeBeginningIndices(index) = appearances(minDistIndex);
    
    appearances(minDistIndex) = [];
    maxMat(1,appearances) = 0;
end
ridgeBeginningIndices = ridgeBeginningIndices(1:index);

%Ridge index variable
rI = 1;

%For each ridge beginning on the lowest scale, follow the ridge and set
%it's properties
for k = ridgeBeginningIndices
    %Get ridge data

    %Set initial length and indices
    ridges(rI).length = 1;
    ridges(rI).initialIndex = k;
    ridges(rI).nextIndex = maxMat(1,k);
    
    for a = 2:numel(SCALES)
        currentInd = ridges(rI).nextIndex;
        
        % Set new maximal CWT coefficient
        currentCoef = abs(coefs(a,currentInd));
        if ridges(rI).maxCoef < currentCoef
            ridges(rI).maxCoef = currentCoef;
            ridges(rI).maxCoefIndex = currentInd;
            ridges(rI).maxCoefScale = a;
        end
        
        nextInd = maxMat(a, currentInd);
        
        if nextInd < 0
           break; 
        elseif nextInd >= currentInd - a && nextInd <= currentInd + a
            %Set new index and length. Then proceed.
            ridges(rI).nextIndex = nextInd;
            ridges(rI).length = ridges(rI).length + 1;
            
            %This position in maxMat has been visited, so it is shut down
            maxMat(a,currentInd) = -1;
        else
            %Continue with next ridge
            break;
        end
    end
    rI = rI + 1;
end

%% Define SNR

%Compute Noise values

%Set lower and upper bound for noise calculation
window = NOISE_WINDOW_SIZE;
ridgeCount = rI;
for r = 1:ridgeCount
    if strcmpi(METHOD,'maxCoef')
        index = ridges(r).maxCoefIndex;
    else
        index = ridges(r).initialIndex;
    end
    
    %Compute window size
    if index > window
        lowerBound = index - window;
    else
        lowerBound = 1;
    end
    if index < numel(spectrum) - window
        upperBound = index + window;
    else
        upperBound = numel(spectrum);
    end
            
    %Compute noise and supply minimal value, if noise is too small
    noise = prctile(abs(coefs(1,lowerBound:upperBound)),95);
    if noise < NOISE_THRESHOLD
        noise = NOISE_THRESHOLD;
    end
        
    ridges(r).SNR = ridges(r).maxCoef/noise;
end

%% Sort out unwanted peaks and convert to an array
peakVec = zeros(1,numel(spectrum));
sNRVec = zeros(1,numel(spectrum));
ridgeLengthVec = zeros(1,numel(spectrum));
index = 0;

%Check if the ridge fullfills all given conditions and save it's
%location in the output array
for r = 1 : ridgeCount
    MaxAmpCheck = false;
    SNRCheck = false;
    LengthCheck = false;
    if ridges(r).length > 1
        % (1) The scale corresponding to the maximum amplitude on the ridge line
        % should be within a certain range
        if ridges(r).maxCoefScale >= MIN_AMP_SCALE && ridges(r).maxCoefScale <= MAX_AMP_SCALE
            MaxAmpCheck = true;
        end
        
        % (2) SNR should be bigger than a certain threshold
        if ridges(r).SNR > SNR_THRESHOLD
            SNRCheck = true;
        end
        % (3) The length of the ridge line should be larger than a certain
        % threshold
        if ridges(r).length > LENGTH_THRESHOLD
            LengthCheck = true;
        end
        
         if MaxAmpCheck && SNRCheck && LengthCheck
             index = index + 1;
             %Save this ridge
             if strcmpi(METHOD,'maxCoef')
                 peakVec(index) = ridges(r).maxCoefIndex;
             else
                 peakVec(index) = ridges(r).initialIndex;
             end
             
             sNRVec(index) = ridges(r).SNR;
             ridgeLengthVec(index) = ridges(r).length;
             
         end
    end
end

%Extract nonzero elements
peakVec = peakVec(1:index);
sNRVec= sNRVec(1:index);
ridgeLengthVec = ridgeLengthVec(1:index);

end

