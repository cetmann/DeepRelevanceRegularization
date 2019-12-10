classdef MSFeatureData < MSData
  % Store feature vectors for mass spec datasets
  %
  % This class stores feature vectors for mass spec datasets. It is derived
  % from the base class MSData and inherits properties for storing data
  % items (feature vectors) and optional related information (annotations, 
  % positions).
  %
  % Properties (in addition to superclass):
  %   creator: Name of the feature map used for creating the feature data
  %   validItems: Logical vector indicating items with valid feature data
  %   featureInfo: Struct for arbitrary feature information
  %
  % Methods:
  %   MSFeatureData: Constructor
  %   reduceFeatures: Reduces number of features
  %
  % MSFeatureData uses the handle semantic, i.e. when assigning an object
  % of this class to a variable, only a reference to the original object is
  % copied. Use the copy method to create a (shallow) copy.
  
  properties (SetAccess = immutable)
    creator  % Name of the feature map used for creating the feature data
  end
  properties
    % Logical vector indicating which items have been assigned valid feature data 
    validItems = [];      
    % Struct for arbitrary feature information
    featureInfo = struct; 
  end
  
  methods
    function obj = MSFeatureData (data, creator, sourceData, validItems)
      % Constructor
      % obj = MSFeatureData(data, creator, sourceData)
      %   data: Feature data matrix containing feature vectors as rows
      %   creator (optional): Feature map used to create the feature data object
      %   sourceData (optional): Source object, derived from MSData. If 
      %     specified, positions are copied from sourceData.
      %   validItems (optional): Logical vector indicating for which items 
      %     valid feature data has been computed
      
      narginchk(1,4);
      if nargin < 2
        creator = '';
      end
      if nargin < 3
        sourceData = [];
      end
      if nargin < 4
        validItems = [];
      end
      % Check input arguments
      if ~ischar(creator)
        error('creator must be either empty or a string');
      elseif ~isempty(sourceData) && ~isa(sourceData, 'MSData')
        error('sourceData must be either empty or an MSData object');
      elseif ~isempty(validItems) && ...
             ~(isvector(validItems) && (islogical(validItems) || isnumeric(validItems)))
        error('validItems must be a logical or numerical vector');
      end
      
      % Call superclass constructor
      obj@MSData(data);
      % Set creator string
      obj.creator = creator;
      % If source data is given, copy positions. (Don't copy annotations, 
      % as these are typically not used on derived data objects.)
      if ~isempty(sourceData)
        if ~isempty(sourceData.positions)
          obj.setPositions(sourceData.positions);
        end
      end
      % Set validItems vector
      if isempty(validItems)
        obj.validItems = true(obj.numItems, 1);
      else
        obj.validItems = logical(validItems(:));
      end
      % Check consistency
      obj.assert;
    end
    
    function fd = reduceFeatures(obj, numFeatures)
        % Reduces number of features
        % INPUT
%         %   numFeatures: number of resulting features after reduction
%         if numFeatures > obj.dataLength
%             warning('input is larger than actual number of features. Nothing is reduced')
%         end
        fd = MSFeatureData(obj.data(:,1:min(numFeatures, obj.dataLength)), obj.creator, obj, obj.validItems);
    end
    function assert (obj)
      % Check object consistency, abort if inconsistent
      
      % Call superclass assert
      assert@MSData(obj);
      % Check validItems vector
      assert(length(obj.validItems) == obj.numItems, ...
             ['Assertion failed: Length of MSFeatureData.validItems (%d) ', ...
              'must match number of data items (%d)'], ...
             length(obj.validItems), obj.numItems);
    end
  end  
end

