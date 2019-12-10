classdef MSDummyTrigger < MSPreprocessingTrigger
    %Dummy preprocessing class
    methods
        %This method does not produce any changes in the maldi object
        function apply(~, ~)
        end
    end
end