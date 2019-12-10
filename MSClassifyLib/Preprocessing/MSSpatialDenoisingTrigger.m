classdef MSSpatialDenoisingTrigger < MSPreprocessingTrigger
    % Class for the spatial denoising preprocessing of data. The logic of the
    % preprocessing aplication is implemented in the method <apply>.
    %
    % Properties
    %   denMethod: Denoising method
    %   denParam: Denosing parameter
    % Methods
    %   MSSpatialDenoising: constructor
    %   apply: modifies the input data by performing a spacial denosing
    %
    % MSSpatialDenoisingTrigger uses the handle semantic, i.e. when assigning an object
    % of this class to a variable, only a reference to the original object is
    % copied. Use the copy method to create a deep copy.
    properties (SetAccess = immutable)
        denMethod; %Denoising method
        denParam; %Denoising parameter
    end
    methods      
        function obj = MSSpatialDenoisingTrigger(denMethod, denParam)
            % Constructor:
            % obj=MSSpatialDenoising. Creates an object where spatial 
            % denosing is applied with default parameters defined in 
            % MSSpatialDenoising function

            %Input validation
            narginchk(0,2);
            if nargin < 2
                obj.denParam = []; % Default parameters will be used
            else                
                if ~isempty(obj.denParam) && isnumeric(obj.denParam) && ~isscalar(obj.denParam)
                    error('denParam must be a numeric scalar')
                end  
                obj.denParam=denParam;
            end
            if nargin < 1
                obj.denMethod = []; % Default denoising method will be used
            else
                if ~isempty(obj.denMethod)&& ~ischar(obj.denMethod)
                    error('denMethod should be a string indicating a valid method')
                end           
                obj.denMethod=denMethod;
            end
        end
        function apply(obj, maldiData)
            % obj.apply(maldiData)
            % Applies the spatial denoising to maldiData
            MSSpatialDenoising(maldiData, obj.denMethod, obj.denParam);
        end        
    end
end