classdef MSReduceTrigger < MSPreprocessingTrigger
    % Preprocessing trigger reducing data to m/z range or item mask
    % 
    % Properties
    %   mzRange: Output m/z range
    %   itemMask: Logical item mask vector. May be specified in constructor 
    %   or using setItemMask().
    %
    % Methods
    %   MSReduceTrigger: Constructor
    %   setItemMask: Specifiy item mask or label
    %   apply: Perform data reduction
    %
    % MSReduceTrigger uses the handle semantic, i.e. when assigning an 
    % object of this class to a variable, only a reference to the original 
    % object is copied. Use the copy method to create a (shallow) copy.
    
    properties (SetAccess = immutable)
      mzRange;
    end
    
    properties (SetAccess = protected)
      itemMask = [];
    end
    
    methods
      function obj = MSReduceTrigger(mzRange, itemMask)
        % Constructor
        % obj = MSReduceTrigger(mzRange, itemMask): Create trigger object
        %   with given m/z range (min, max) and item mask (label data,
        %   logical or numerical vector). Each one of the arguments may be
        %   omitted, but not both.

        narginchk(1,2);
        obj.mzRange = [];
        obj.itemMask = [];
        isValidRange = @(x) isnumeric(x) && isvector(x) && ...
                            length(x) == 2 && x(2) >= x(1);
        isValidMask = @(x) (isa(x, 'MSLabelData') && x.dataLength == 1) || ...
                           (isvector(x) && (isnumeric(x) || islogical(x)));
        if nargin == 2
          % Both arguments given
          if isempty(mzRange) || isValidRange(mzRange)
            obj.mzRange = mzRange;
          else
            error('Invalid m/z range specified');
          end
          if isempty(itemMask) || isValidMask(itemMask)
            obj.setItemMask(itemMask);
          else
            error('Invalid item mask specified');
          end
        else
          % Only one argument given
          if isempty(mzRange) || isValidRange(mzRange)
            obj.mzRange = mzRange;
          elseif isValidMask(mzRange)
            obj.setItemMask(mzRange);
          else
            error('Invalid argument specified');
          end
        end
      end
      
      function setItemMask (obj, M)
        % Specifiy item mask or label
        % obj.setItemMask(M): Specify item mask M as a logical or numerical 
        %   vector, or as an MSLabelData object.
        
        if islogical(M) || isnumeric(M)
          if ~(isempty(M) || isvector(M))
            error('Item mask must be a logical or numerical vector');
          else
            if isnumeric(M)
              M = M > 0;
            end
            obj.itemMask = M(:)';
          end
        elseif isa(M, 'MSLabelData')
          if M.dataLength > 1
            error('Item mask must be a single column label object');
          else
            obj.itemMask = M.data > 0;
          end
        end
      end

      function apply (obj, maldiData)
        % Apply m/z range and/or data item reduction to maldiData
        
        if ~isempty(obj.itemMask)
          % Reduce to data item subset
          maldiData.reduce(obj.itemMask);
        end
        if ~isempty(obj.mzRange)
          % Reduce m/z range, limit to intersection with data range
          newMzRange = [max(obj.mzRange(1), maldiData.mzVector(1)), ...
                     min(obj.mzRange(2), maldiData.mzVector(end))];
          if newMzRange(2) <= newMzRange(1)
            error('Data m/z range outside specified output range')
          end
          maldiData.reduceMzRange(newMzRange);
        end
      end
    end
end