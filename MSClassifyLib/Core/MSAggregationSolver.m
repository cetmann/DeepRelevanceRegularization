classdef MSAggregationSolver < matlab.mixin.Copyable
    % Abstract class which defines the logic of the aggregation
    % Methods
    %   resolveAggregation: Defines how to aggregate according to the
    %   number of spots for each class in the segment
    %   resolveAggregation_impl: (Abstract)
    
    methods
        function class=resolveAggregation(obj, occurencies, itemsTotal)
        %Defines how to aggregate according to the number of spots for each
        %class in the segment
        % INPUT
        %   occurencies: vector with as many elements as classes indicating
        %               the occurrencies of spots at aggregation level.
        %   itemsTotal: total number of spots for a specific aggregation
        %               level label
        % OUTPUT
        %   class: resulting class of the label
        
        %check input
        if ~isvector(occurencies)||~isnumeric(occurencies)||any(occurencies<0)
            error('The input <occurencies> must be a numerical non-negative vector')
        end
        if ~isnumeric(itemsTotal)||~isscalar(itemsTotal)||itemsTotal<sum(occurencies)
            error('The input <itemsTotal> must be a positive number')
        end
        occurencies=floor(occurencies);
        itemsTotal=floor(itemsTotal);
        %call specific implementation of aggregation
        class=obj.resolveAggregation_impl(occurencies, itemsTotal);
        end
    end
    methods (Abstract, Access=protected)
        % aggregation logic
        class=resolveAggregation_impl(obj, occurencies, itemsTotal)
    end
end