classdef MSMeanFilteringTrigger < MSPreprocessingTrigger
    % Class for the mean filtering of the data. The logic of the
    % filtering aplication is implemented in the method <apply>
    %
    % Methods
    %   apply: Contains the filtering logic
    %
    % MSMeanFilteringTrigger uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    methods         
        function apply(~, maldiData)
            %obj.apply(maldiData)
            %smooth data using mean filtering implemented in MSSmooth
            MSSmooth(maldiData);
        end        
    end
end