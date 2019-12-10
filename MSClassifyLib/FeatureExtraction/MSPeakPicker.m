classdef MSPeakPicker < matlab.mixin.Copyable
    % Peak picking based on the MATLAB function mspeaks.
    %
    % This class implements a peak picking method based on the MATLAB function
    % mspeaks. In a given spectrum, the n highest peaks satisfying certain criteria are
    % sought.
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
    %    MSPeakPicker: Constructor
    %    createPeakList: Create a list containing the n most important peaks
    %    show: Show a plot in which the peaks from peakList are highlighted
    %    meanDist: Calculate mean distance of all peaks in obj.spectrum
    %
    % MSPeakPicker uses the handle semantic, i.e. when assigning an object
    % of this class to a variable, only a reference to the original object is
    % copied. Use the copy method to create a deep copy.
    
    properties (SetAccess = protected)
        spectrum
        maldiData
        peakList
        peakIndices
    end
    
    methods
        function obj = MSPeakPicker(pMaldiData, pSpectrum)
            % Constructor
            % obj = MSPeakPicker(pMaldiData, spec):
            %  pMaldiData: MSData object containing the m/z vector belonging to
            %  the spectrum
            %  pSpectrum: Spectrum in which the peaks will be picked
            
            %Check input arguments and initialize properties
            if ~isnumeric(pSpectrum) || isempty(pSpectrum)
                throw(MException('MSClassifyLib:MSPeakPicker:EmptyOrNonNumericSpectrumInput',...
                    ['The second input variable for MSPeakPicker() must be non-empty and of type numeric.']));
            end
            obj.spectrum = pSpectrum;
            
            if ~isa(pMaldiData, 'MSData')
                throw(MException('MSClassifyLib:MSPeakPicker:InvalidMaldiDataInput',...
                    'The first input variable for MSPeakPicker() must be an MSData object'));
            end
            obj.maldiData = pMaldiData;
            
            obj.peakList = [];
            obj.peakIndices = [];
        end
        
        function output =  createPeakList(obj, n, varargin)
            % This method extracts the n most important peaks from the spectrum.
            %
            % --- Input parameters ---
            % integer n: Number of peak indices to find (Mandatory)
            % 'HeightFilter',float hf: The minimum height for a peak in the spectrum to be
            % selected. Please enter values in the units of the m/z vector.(Optional)
            % 'OverSegmentation',float os: Minimum distance between neighboring peaks.
            % Please enter values in the units of the m/z vector. (Optional)
            %
            % --- Output ---
            % peakList: A list containing the m/z values at which the peaks are
            % located
            
            
            %Create local variables with shorter names
            spec = obj.spectrum;
            mz = obj.maldiData.mzVector;
            
            %Check input for n
            if n < 1 || ~isnumeric(n) || numel(n) > 1 || isempty(n) ||...
                    mod(n,1) ~= 0  || n > numel(spec)
                error('The input variable n must be a positive integer, smaller than numel(spectrum).');
            end
            
            % Check if spectrum and mzs have the same size – this is needed for
            % mspeaks() later.
            if size(spec) ~= size(mz)
                mz = mz';
            end
            
            % Declare default values
            % The real default values are calculated late, so that calculation
            % only happens, if the default values are used.
            overSegmentation = [];
            heightFilter = [];
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
                        case 'oversegmentation'
                            input = varargin{iArgument + 1 };
                            if isnumeric(input) && input > 0 && numel(input) == 1
                                overSegmentation = input;
                            else
                                warning('The input for OverSegmentation is invalid. The default value 10*meanDist()+1 will be used.');
                            end
                            iArgument = iArgument + 2;
                            
                        case 'heightfilter'
                            input = varargin{iArgument + 1 };
                            if isnumeric(input) && input > 0 && numel(input) == 1
                                heightFilter = input;
                            else
                                warning('The input for HeightFilter is invalid. The input will be ignored.');
                            end
                            iArgument = iArgument + 2;
                            
                        otherwise
                            %Handle incorrect input
                            iArgument = iArgument + 1;
                            error(['Unknown argument! The argument at input position ' num2str(iArgument) ' is invalid.']);
                    end
                end
            end
            if isempty(overSegmentation)
                overSegmentation = 10 * obj.meanDist() + 1;
            end
            if isempty(heightFilter)
                heightFilter = median(obj.spectrum);
            end
            
            % --- Generate and rank output ---
            output = mspeaks(mz,spec,...
                'OverSegmentationFilter',overSegmentation,...
                'HeightFilter',heightFilter);
            
            % Sort mspeaks output descending in peak height
            [~,sortInd] = sort(output(:,2),'descend');
            output = output(sortInd,1);
            
            %Extract n most important peaks, if possible.
            if n <= numel(output)
                output = output(1:n);
            else
                warning('The count of peaks that were found, is less than n.');
            end
            
            output = output';
            
            %For every value in output, find the matching index and save it in
            %peakIndices
            obj.peakIndices = zeros([1,numel(output)]);
            for k = 1:length(output)
                [~,ind] = min(abs(mz-output(k)));
                obj.peakIndices(k) = ind;
            end
            obj.peakList = output;
        end
        
        function figureHandle =  show(obj,varargin)
            % Show a plot of the spectrum, in which the peaks from the peak list
            % are marked.
            %
            % -- Input parameters ---
            % 'Title',string title - A string containing an optional title for the plot
            %
            % --- Output ---
            % f - A handle to the figure conaining the plot
            
            if isempty(obj.peakList) || isempty(obj.peakIndices)
                error('Please create a peakList before calling this objects show method.');
            else
                mz = obj.maldiData.mzVector;
                spec = obj.spectrum;
                
                %Create figure spanning the whole screen
                figureHandle = figure('Units','Normalize','Position',[0 0 1 1]);
                hold on
                
                %Plot spectrum
                plot(mz,spec);
                axis tight;
                
                %Plot peak markers
                plot(mz(obj.peakIndices),spec(obj.peakIndices),'ro','MarkerSize',8,'LineWidth',1.4);
                
                %Set axis properties
                ax = gca;
                ax.XColor = repmat(.55,[1,3]);
                ax.YColor = repmat(.55,[1,3]);
                ax.Box = 'off';
                ax.TickLength = [0 0];
                ax.FontSize = 10;
                
                %Set title if necessary
                if nargin > 1
                    title(varargin{1});
                end
                
            end
        end
        
        function mD = meanDist(obj)
            %Calculate mean distance of two peaks in the spectrum
            [~,localMaxIndices] = findpeaks(double(obj.spectrum));
            localMaxIndices = obj.maldiData.mzVector(localMaxIndices);
            mD = ceil(mean(diff(localMaxIndices)));
        end
        
        function p = get.peakList(obj)
            p = obj.peakList;
        end
        
        function s = get.spectrum(obj)
            s = obj.spectrum;
        end
        
        function obj = setSpectrum(obj,spec)
            if isempty(spec) || ~isnumeric(spec) || numel(spec) ~= obj.maldiData.dataLength
                throw(MException('MSClassifyLib:MSPeakPicker:EmptyOrNonNumericSpectrumInputForSetterMethod',...
                    ['The input spectrum must be a non-empty, numeric array containing ' num2str(numel(obj.maldiData.mzVector)) ' items.']));
            else
                obj.spectrum = spec;
            end
        end
    end
end

