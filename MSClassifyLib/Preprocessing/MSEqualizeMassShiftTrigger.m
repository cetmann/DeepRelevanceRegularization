classdef MSEqualizeMassShiftTrigger < MSPreprocessingTrigger
    % Preprocessing trigger performing mass shift equalization
    % 
    % Properties
    %   args:  Parameters (see MSEqualizeMassShift)
    %   prfIn: Reference profile
    %
    % Methods
    %   MSEqualizeMassShiftTrigger: Constructor
    %   setReferenceProfile: Set reference profile for MSEqualizeMassShift
    %   apply: Perform mass shift equalization
    %
    % MSEqualizeMassShiftTrigger uses the handle semantic, i.e. when
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a
    % (shallow) copy.
    
    properties (SetAccess = immutable)
      args;  % Parameters (see MSEqualizeMassShift)
      prfIn; % Reference profile
    end
    
    methods
      function obj = MSEqualizeMassShiftTrigger (varargin)
        % Constructor
        % obj = MSEqualizeMassShiftTrigger(): Create trigger object using
        %   default parameters
        % obj = MSEqualizeMassShiftTrigger(arg1, arg2, ...): Specify
        %   additional arguments to be passed to MSEqualizeMassShift()

        % Store input arguments
        obj.args = varargin;
        obj.prfIn = [];
      end
      
      function setReferenceProfile (obj, P)
        % Set reference profile
        % obj.setReferenceProfile(P): Specify reference profile P as
        %   supported by MSEqualizeMassShift(). 
        % Note: This will not work if a reference profile has already been
        % specified in the constructor.
        
        obj.prfIn = P;
      end

      function apply (obj, maldiData)
        % If reference profile is specified, combine with argument list
        A = obj.args;
        if ~isempty(obj.prfIn)
          A = [{obj.prfIn}; A(:)];
        end
        % Apply mass shift equalization to maldiData
        msdOut = MSEqualizeMassShift(maldiData, A{:});
        % Store transformed MALDI data in input object
        maldiData.initFrom(msdOut);
      end
    end
end