classdef MSSqrtTransformationTrigger < MSPreprocessingTrigger
    % Class for the sqrt transformation of the data. The logic of the
    % preprocessing aplication is implemented in the method <apply>
    % 
    % Properties
    % 
    % Methods
    %   MSSqrtTransformationTrigger: Constructor
    %   apply: Contains logic of the transformation
    %
    % MSLogTransformationTrigger uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    
    properties (SetAccess = immutable)
        nPeak; 
        nPeakLogIntensitySlope; 
    end
    
    methods        
        function obj = MSSqrtTransformationTrigger()
        end
        
        function apply(obj, maldiData) %#ok<INUSL>
            % obj.apply(maldiData)
            % applies the transformation to the MSMaldiData maldiData
            sqrt_transformation( maldiData );
        end
        
    end
end

function [ data ] = sqrt_transformation( D )
    nargoutchk(0,1);
    if nargout > 0
      data = D.copy;
      data.data = sqrt(1+max(0,data.data));
    else
      D.data = sqrt(1+max(0,D.data));
    end
end
