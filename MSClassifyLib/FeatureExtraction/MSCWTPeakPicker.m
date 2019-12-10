classdef MSCWTPeakPicker < MSPeakPicker
    % Peak picking based on the function cWTPeaks. 
    %
    % This class implements a peak picking method based on the function
    % cWTPeaks. In a given spectrum, the n highest peaks satisfying certain criteria are
    % sought. It is a subclass of MSPeakPicker and inherits all of it's
    % properties and methods. The only difference is, that this class has a
    % different 'createPeakList' method.
    %
    % Properties (read-only):
    %    maldiData: A MSData object holding the m/z vector belonging to the spectrum
    %    peakList: A list of peaks selected by the the method createPeakList
    %    peakIndices: The indices belonging to the peaks in peakList
    %
    %    (writeable):
    %    spectrum: The spectrum in which peaks will be sought
    %
    % Methods:
    %    MSCWTPeakPicker: Constructor
    %    createPeakList: Create a list containing the n most important peaks
    %    show: Show a plot in which the peaks from peakList are highlighted
    %    meanDiff: Calculate mean difference of all peaks in obj.spectrum
    %
    % MSCWTPeakPicker uses the handle semantic, i.e. when assigning an object
    % of this class to a variable, only a reference to the original object is
    % copied. Use the copy method to create a deep copy.
     
     
     methods
         function obj = MSCWTPeakPicker(pMaldiData, pSpectrum)
             % Constructor
             % obj = MSCWTPeakPicker(pMaldiData, spec):
             %  pMaldiData: MSData object containing the m/z vector belonging to
             %  the spectrum
             %  pSpectrum: Spectrum in which the peaks will be picked
             if nargin == 0
                 pMaldiData = [];
                 pSpectrum = [];
             end
             obj@MSPeakPicker(pMaldiData,pSpectrum);
         end

        function output =  createPeakList(obj, n, varargin)
           
            % This method extracts the n most important peaks from the spectrum using the function cWTPeaks.
            % 
            % For more information please see the documentation of cWTPeaks
            %
            % --- Input parameters ---
            % n: Number of peak indices to find (Mandatory)
            % Wavelet: The type of wavelet to use
            % lengthThreshold: The minimum ridge lenght
            % SNRThreshold: The minimum Signal-To-Noise-Ratio value
            % 
            
            
            %Create abbreviations to improve code readability
            spec = obj.spectrum;
            
            %Check input for n
            if n < 1 || ~isnumeric(n) || numel(n) > 1 || isempty(n) ||...
                    mod(n,1) ~= 0  || n > numel(spec)
                error('The input variable n must be a positive integer, smaller than numel(spectrum).');
            end
                       
          
            %Initialze additional input variables
            wavelet = 'mexh';
            lengthThreshold = 1;
            sNRThreshold = .1;
            scales = 1:12;
            
            % --- Check and process additional input variables ---
            mandatoryArguments = 2;
            if nargin > mandatoryArguments
                iArgument = 1;
                while iArgument <= nargin - mandatoryArguments
               
                    currArg = varargin{iArgument};
                    if isempty(currArg) || isnumeric(currArg)
                        currArg = '';
                    end
                    
                    switch lower(currArg)
                        case 'wavelet'
                            input = varargin{iArgument + 1};
                            if ~ischar(input)
                               warning('The input for Wavelet is invalid. It must be of type char. A default value will be used.');
                            else
                               wavelet = input; 
                            end
                            
                            iArgument = iArgument + 2;
                    
                        case 'lengththreshold'
                            input = varargin{iArgument + 1};
                            if isnumeric(input) && input > 0 && numel(input) == 1
                                lengthThreshold = input;
                            else
                                warning('The input for LenghtThreshold is invalid. It must be a positive integer. A default value will be used.');
                            end
                            iArgument = iArgument + 2;
                        
                        case 'snrthreshold'
                            input = varargin{iArgument + 1};
                            if isnumeric(input) && input > 0 && numel(input) == 1
                                sNRThreshold = input;
                            else
                                warning('The input for SNRThreshold is invalid. It must be a positive integer. A default value will be used.');
                            end
                            iArgument = iArgument + 2;
                            
                        case 'scales'
                            input = varargin{iArgument + 1};
                            if isnumeric(input) && min(size(input)) == 1
                                scales = input;
                            else
                                warning('The input for Scales is invalid. It must be an integer array. A default value will be used.');
                            end
                            iArgument = iArgument + 2;
                            
                        otherwise
                            %Handle incorrect input
                            iArgument = iArgument + 1;
                            error(['Unknown argument! The argument at input position ' num2str(iArgument) ' is invalid.']);
                    end
                end
            end
            
            % --- Generate and rank output ---
            [output, SNR, ridgeLength] =  cWTPeaks(spec,...
                'Wavelet',wavelet,...
                'LengthThreshold',lengthThreshold,...
                'SNRThreshold',sNRThreshold,...
                'Scales',scales);
            
            rank = SNR .* ridgeLength;
            [~,sortInd] = sort(rank,'descend');
            output = output(sortInd);
                        
          
            %Extract n most important peaks, if possible.
            if n <= numel(output)
                output = output(1:n);
            else
                warning(['The count of peaks that were found, is less than ' num2str(n) '.']);
            end
            
            %Set properties
            obj.peakIndices = output;
            output = obj.maldiData.mzVector(output);
            obj.peakList = output;
        end
    end
end

