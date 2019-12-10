classdef MSLogTransformationTrigger < MSPreprocessingTrigger
    % Class for the log transformation of the data. The logic of the
    % preprocessing aplication is implemented in the method <apply>
    % 
    % Properties
    %   nPeak:     The nPeak highest peak of the max-spectrum of the data
    %     will be transformed with the log transformation 
    %       f(x) = b*log(1+x/b) 
    %     such that
    %       f'(x_nPeak) = nPeakLogIntensitySlope = p.
    %     b will be chosen accordingly, i.e. b := p/(1-p)*x_nPeak.
    %     The above function fulfills f'(0) = 1 for all b.
    %   nPeakLogIntensitySlope: compare above
    % Methods
    %   MSLogTransformationTrigger: Constructor
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
        function obj = MSLogTransformationTrigger(nPeak,nPeakLogIntensitySlope)
          % Constructor  
          % obj = MSLogTransformationTrigger()
          % returns an instance with the specified normalization logic
          
          narginchk(0,2);
          if nargin < 1 || isempty(nPeak)
            nPeak = 30;
          end
          if nargin < 2 || isempty(nPeakLogIntensitySlope)
            nPeakLogIntensitySlope = 0.01;
          end
          obj.nPeak = nPeak;
          obj.nPeakLogIntensitySlope = nPeakLogIntensitySlope;
        end
        
        function apply(obj, maldiData)
            % obj.apply(maldiData)
            % applies the transformation to the MSMaldiData maldiData
            log_transformation( maldiData, obj.nPeak, obj.nPeakLogIntensitySlope );
        end
        
    end
end

function [ data ] = log_transformation( D, nPeak, nPeakLogIntensitySlope )
    narginchk(1,3);
    nargoutchk(0,1);
    if nargin < 3
      nPeakLogIntensitySlope = 0.05; % 0.01
    end
    if nargin < 2
      nPeak = 30;
    end
    
    p = nPeakLogIntensitySlope; % short notation
    mx = max(0,max(D.data,[],1));
    [~,I] = sort(mx,'descend'); %mx(I(1)) = max(mx)
    revI = 1:length(I);
    revI(I) = 1:length(I); %I(revI(j)) = j, that is we can delete j from I by I(revI(j)) = []
    mask = false( 1, length(I) );
    for i=1:nPeak
        mxIdx = I( find(~mask,1) );
        mask(revI((D.mzVector(mxIdx)-0.5 <= D.mzVector) & (D.mzVector <= D.mzVector(mxIdx)+0.5))) = 1;
    end
    b = p/(1-p)*mx(I(1));
    
    if nargout > 0
      data = D.copy;
      %data.data = b*log( 1+max(0,data.data)./ b );
      data.data = b*log1p( max(0,data.data)./ b ); % Numerically robust for small x
    else
      %D.data = b*log( 1+max(0,D.data)./b );
      D.data = b*log1p( max(0,D.data)./ b ); % Numerically robust for small x
    end
end
