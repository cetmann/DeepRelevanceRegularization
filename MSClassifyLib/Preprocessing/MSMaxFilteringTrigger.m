classdef MSMaxFilteringTrigger < MSPreprocessingTrigger
    % Class for the max filtering of the data. The logic of the
    % filtering aplication is implemented in the method <apply>
    % 
    % Properties
    %   se: structural element used for filtering (see function imdilate)
    % Methods
    %   MSMaxFilteringTrigger: Constructor
    %   apply: Contains logic of the filtering
    %
    % MSMaxFilteringTrigger uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    properties (SetAccess = immutable)
        se; %structural element used for filtering (see function imdilate)
    end
    methods   
        function obj = MSMaxFilteringTrigger(se)
            %Constructor
            %obj.MSMaxFilteringTrigger-uses the default structural element
            %[1 1 1] for performing filtering
            %obj.MSMaxFilteringTrigger(se)-uses the structural element se
            %for performing filtering
           
            narginchk(0,1)
            if nargin < 1
                obj.se=[1 1 1];
            end
            obj.se=se;
        end
        function apply(obj, maldiData)
            %obj.apply(maldiData)
            % applies the max filtering to maldiData implemented in imdilate
            maldiData.data = imdilate(maldiData.data, obj.se);
        end        
    end
end