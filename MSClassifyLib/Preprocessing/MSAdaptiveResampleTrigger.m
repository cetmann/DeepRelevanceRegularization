classdef MSAdaptiveResampleTrigger < MSPreprocessingTrigger
    % Preprocessing trigger performing adaptive resampling
    % 
    % Properties
    %   params:   Parameter struct (see MSAdaptiveResampleParams)
    %   equalize: If true, use mass shift equalization on single spectra
    %   itemMask: Logical item mask matrix defining data item subsets for
    %             which adaptive resampling is performed indepedently. May
    %             be specified in constructor or using setItemMask(), not
    %             allowed with 'equalize'.
    %
    % Methods
    %   MSAdaptiveResampleTrigger: Constructor
    %   setItemMask: Specifiy item mask or label
    %   apply: Perform adaptive resampling
    %
    % MSAdaptiveResampleTrigger uses the handle semantic, i.e. when
    % assigning an object of this class to a variable, only a reference to
    % the original object is copied. Use the copy method to create a
    % (shallow) copy.
    
    properties (SetAccess = immutable)
      params;
    end
    
    properties (SetAccess = protected)
      equalize = false;
      itemMask = [];
    end
    
    methods
      function obj = MSAdaptiveResampleTrigger(varargin)
        % Constructor
        % obj = MSAdaptiveResampleTrigger(): Create trigger object using
        %   default parameters (see MSAdaptiveResampleParams)
        % obj = MSAdaptiveResampleTrigger('equalize'): Use mass shift
        %  equalization to process each spectrum separately
        % obj = MSAdaptiveResampleTrigger(itemMask): Specify item subsets
        %   on which adaptive resampling is performed independently
        % obj = MSAdaptiveResampleTrigger(__, name, value, ...): Specify
        %   parameters as name value pairs (see MSAdaptiveResampleParams)
        
        if nargin >= 1 
          if ~ischar(varargin{1})
            % First argument is not a parameter name, must be itemMask
            obj.setItemMask(varargin{1});
            varargin(1) = [];
          elseif strcmpi(varargin{1}, 'equalize')
            obj.equalize = true;
            varargin(1) = [];
          end
        end
        obj.params = MSAdaptiveResampleParams(varargin{:});
      end
      
      function setItemMask (obj, M)
        % Specifiy item mask or label
        % obj.setItemMask(M): Specify item mask M as a logical vector or
        %   matrix, as a numerical vector, or as an MSLabelData object.
        %   If specified, adaptive resampling is performed separately for
        %   each item subset.
        
        if obj.equalize && ~isempty(M)
          error('Item mask must not be specified with ''equalize''');
        end
        obj.itemMask = MSCheckItemMask([], M);
      end

      function apply (obj, maldiData)
        % Apply adaptive resampling to maldiData
        args = obj.paramsToArgs();
        if obj.equalize
          args = [{'equalize'}; args];
        else
          args = [{obj.itemMask}; args];
        end
        msdOut = MSAdaptiveResample(maldiData, args{:});
        % Store resampled MALDI data in input object
        maldiData.initFrom(msdOut);
      end
    end
    
    methods (Access = protected)
      function args = paramsToArgs (obj)
        % Convert parameter struct to name-value argument list
        names = fieldnames(obj.params);
        values = struct2cell(obj.params);
        args = cell(2*length(names), 1);
        args(1:2:end-1) = names;
        args(2:2:end) = values;
      end
    end
end