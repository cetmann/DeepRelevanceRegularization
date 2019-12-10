classdef MSNormalizationTrigger < MSPreprocessingTrigger
    % Class for the normalization of the data. The logic of the
    % preprocessing aplication is implemented in the method <apply>
    % 
    % Properties
    %   normalization: String indicating the normalization type
    %   init:          Init flag
    % Methods
    %   MSNormalizationTrigger: Constructor
    %   apply: Contains logic of the normalization
    %
    % MSNormalizationTrigger uses the handle semantic, i.e. when 
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a deep copy.
    
    properties (SetAccess = immutable)
        normalization; % String indicating the normalization type
        init;          % Init flag. If this set to true, the trigger 
                       % re-computes the normalization, otherwise the
                       % existing normalization is restored (the default)
    end
    methods        
        function obj = MSNormalizationTrigger(normalization, init)
          % Constructor  
          % obj=MSNormalizationTrigger(normalization)
          % returns an instance with the specified normalization logic
          
          %input validation
          narginchk(0,2)
          if nargin<1
              %default value
              obj.normalization='tic';
          else
              if ~ischar(normalization) 
                  error('normalization must be a string')
              end
              obj.normalization=normalization;
          end
          if nargin<2
              obj.init=false;
          else
              if ~(islogical(init) && isscalar(init))
                  error('init must be a logical scalar');
              end
              obj.init=init;
          end
          %checking whether the normalization string is a valid
          %normalization is done in setNormalization in the apply
          %method so that the code is more flexible to future inclusions 
          %of normalization types 
        end
        
        function apply(obj, maldiData)
            % obj.apply(maldiData)
            % applies the normalization to the MSMaldiData maldiData
            if obj.init
              maldiData.initNormalization;
            end
            maldiData.setNormalization(obj.normalization)
        end
        
    end
end