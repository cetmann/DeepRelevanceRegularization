classdef MSWeightedAggregationSolver < MSAggregationSolver
    %This class defines the logic for resolving the aggregation using a
    %weight vector
    %
    % Properties
    %   weights: vector of weights used for aggregation
    %   
    % Methods
    %   MSWeightedAggregationSolver: Constructor
    %   resolveAggregation: returns a label (scalar) to a specific 
    
    properties (SetAccess=immutable)
        weights; %vector of weights used for aggregation
    end
    
    methods
        function obj=MSWeightedAggregationSolver(weights)
            %Constructor
            % INPUT
            %   weights: numerical vector with weights for the different
            %   classes, which are applied in the aggregation process
            if ~isvector(weights)||~isnumeric(weights)
                error('The input argument must be a numeric vector')
            end
            obj.weights=weights(:);
        end
    end
    methods (Access=protected)
        function class=resolveAggregation_impl(obj, occurencies, ~)
        % uses the weights to chose the weighted maximum
        if length(occurencies)~=length(obj.weights)
            error(['The number of classes (length of occurencies) must'...
                 'coincide with the length of the property <weights>'])
        end
        if(occurencies(1)*occurencies(2)~=0)
            a=1;
        end
        [~, class]=max(occurencies(:).*obj.weights);
        end
    end
end
   