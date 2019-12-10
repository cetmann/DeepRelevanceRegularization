classdef MSVarianceNormalizationTrigger < MSPreprocessingTrigger
    % Class for the variance normalization of data. The logic of the
    % preprocessing aplication is implemented in the method <apply>
    % 
    % Properties
    %   varEstType:     String indicating the var estimation type
    %   init:           Init flag
    % Methods
    %   MSVarianceNormalizationTrigger: Constructor
    %   apply: Contains logic of var normalization
    %
    % MSVarianceNormalizationTrigger uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    
    properties (SetAccess = immutable)
        varEstType; % String indicating the variance estimation type
    end
    properties (Hidden, SetAccess = immutable)
        allowedVarEstTypes = {'mad','std','std0'};
    end
    
    
    methods        
        function obj = MSVarianceNormalizationTrigger( varEstType )
            % Constructor  
            % obj=MSVarianceNormalizationTrigger( varEstType )
            % returns an instance with the specified varEstType logic

            %input validation
            narginchk(0,1)
            if nargin<1
                %default value
                obj.varEstType = 'mad';
            else
                dummyObj = MSVarianceNormalizationTrigger();
                if ~ischar( varEstType ) || ~any(strcmp(varEstType,dummyObj.allowedVarEstTypes))
                    error('Specify varEstType as one of the following: %s', ...
                        strjoin( dummyObj.allowedVarEstTypes, ', ') );
                end
                obj.varEstType = varEstType;
            end
        end
        
        function apply (obj, fData)
        % Apply variance normalization to data
            varEstFun = obj.getVarEstFun( obj.varEstType );
            varEst = varEstFun( fData.data );
            meanEst = mean( fData.data, 1 );
            var0Mask = varEst./abs(meanEst) < eps;
            if any( var0Mask )
                warning(['Certain variables close to 0 variance ',...
                    '(their var normalization has been skipped)']);
            end
            varEst( var0Mask ) = 1;
            M = ones( fData.numItems, 1 )*(1./varEst);            
            fData.data = M.*fData.data;
        end
    end
    
    methods (Static)
        function varEstFun = getVarEstFun( type )
            switch type
                case 'mad'
                    varEstFun = @(X)mad(X,1,1);
                case 'std'
                    varEstFun = @(X)std(X,0,1);
                case 'std0'
                    varEstFun = @(X)sqrt(mean(X.^2,1));
            end
        end
    end
    
end

