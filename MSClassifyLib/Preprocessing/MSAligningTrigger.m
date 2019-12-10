classdef MSAligningTrigger < MSPreprocessingTrigger
    % Class for the alignment of the data. The logic of the
    % alignment is implemented in the method <apply>
    % 
    % Properties
    %   nRefPeaks: Number of reference peaks (integer) 
    %   minMzDist: Minimum mz distance (positive real number)
    %
    % Methods
    %   MSAligningTrigger: Constructor
    %   apply: Contains logic of the alignment
    %
    % MSAligningTrigger uses the handle semantic, i.e. when assigning an object
    % of this class to a variable, only a reference to the original object is
    % copied. Use the copy method to create a deep copy.
    properties (SetAccess = immutable)
        nRefPeaks;
        minMzDist;
    end
    methods
        function obj = MSAligningTrigger(nRefPeaks, minMzDist)
            %Constructor
            %obj.MSAligningTrigger-uses default nRefPeaks and minMzDist
            %obj.MSAligningTrigger(nRefPeaks)
            %obj.MSAligningTrigger(nRefPeaks,[])-use default minMzDist
            %obj.MSAligningTrigger()||obj.MSAligningTrigger([],[])-uses
            %       both default parameters
            
            %input validation
            narginchk(0,2)
            if nargin < 2
                obj.minMzDist=[];
            end
            if nargin < 1
                obj.nRefPeaks=[];
            end
            if ~isempty(minMzDist)&&(~isscalar(minMzDist) || ~isnumeric(minMzDist) || minMzDist<0)
                 error('minMzDist must be a real non-negative number')
            end
            obj.minMzDist=minMzDist;
            if ~isempty(nRefPeaks)&&(~isscalar(nRefPeaks) || ~isnumeric(nRefPeaks) || nRefPeaks<1)
                    error('minMzDist must be a positive integer')
            end         
            obj.nRefPeaks=floor(nRefPeaks);                
        end
        
        function apply(obj, maldiData)
            %obj.apply(maldiData)
            %Applies the logic of peak alignment implemented in MSAlign
            MSAlign(maldiData, obj.nRefPeaks, obj.minMzDist);
        end
    end
end