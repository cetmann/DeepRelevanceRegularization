classdef MSPreprocessingTrigger < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    % Base class for the preprocessing trigger hierarchy. An
    % MSPreprocessingTrigger contains a preprocessing logic implemented in
    % the apply. The apply method is abstract in this class and should be
    % implemented in classes deriving from MSPreprocessingTrigger
    %
    % Methods
    %   apply: (Abstract) Applies the preprocessing logic, modifying maldiData.
    %
    % MSPreprocessingTrigger uses the handle semantic, i.e. when assigning an object
    % of this class to a variable, only a reference to the original object is
    % copied. Use the copy method to create a deep copy.
    
    methods (Abstract)
        % Applies the preprocessing logic, modifying maldiData.
        % For the sake of performance the validation of the input maldiData
        % as MSMaldiData is skipped (it is checked in the method client 
        % of the class)
        apply(obj, maldiData) 
    end
end