classdef MSDataPartition
  % Store a set of labels for a common dataset
  % 
  % This class stores a set of MSLabelData objects associated with a common
  % base dataset object (MSData or MSMaldiData). MSDataPartition objects
  % are used to provide data partition information for a certain
  % classification task.
  %
  % Properties:
  %   The following properties are MSLabelData objects associated with the
  %   common base dataset object. Only classes need to be defined for a
  %   given classification task, other labels are optional.
  %
  %   numItems: Number of items in each label
  %   classes: Ground truth class assignment
  %   samples: Tissue sample / biopsy / core ID
  %   individuals: Individual / patient ID
  %   sampleSets: Sample set / TMA label
  %   cvEntities: Cross validation entities to which cross validation
  %               segmentation must be aligned
  %   other: Struct with any other label data
  %
  %   specificLabels: Cell string array of specific label data properties
  %   otherLabels: Cell string array of label data in 'other' struct
  %
  % Methods:
  %   reduce: Reduce label data to specified item subset
  %   removeUnused: Remove unused labels
  %   assert: Check object consistency, abort if inconsistent
  
  properties
    classes = [];     % Ground truth class assignment
    samples = [];     % Tissue sample / biopsy / core ID
    individuals = []; % Individual / patient ID
    sampleSets = [];  % Sample set / TMA label
    cvEntities = [];  % Cross validation entities to which cross validation
                      % segmentation must be aligned
    other = struct;   % Struct with any other label data
  end

  properties (Dependent)
    numItems;         % Number of items in each label
    specificLabels;   % Cell string array of specific label data properties
    otherLabels;      % Cell string array of label data in 'other' struct
  end

  methods
    function N = get.numItems (obj)
      if isempty(obj.classes)
        N = [];
      else
        N = obj.classes.numItems;
      end
    end
    
    function S = get.specificLabels (obj)
      names = {'classes', 'samples', 'individuals', 'sampleSets', 'cvEntities'};
      S = names(cellfun(@(x) isa(obj.(x), 'MSLabelData'), names));
    end
    
    function S = get.otherLabels (obj)
      names = fieldnames(obj.other);
      S = names(cellfun(@(x) isa(obj.other.(x), 'MSLabelData'), names));
    end
    
    function P = applyLabelFunction (obj, F)
      % Apply specified function to all labels
      % P = obj.applyLabelFunction(F): Apply function F to all labels in
      %   obj and store resulting labels in new MSDataPartition P.

      % Create new output object
      P = MSDataPartition();
      % Call F on defined properties
      properties = obj.specificLabels;
      for k = 1:length(properties)
        P.(properties{k}) = F(obj.(properties{k}));
      end
      % Call F on other label data members
      properties = obj.otherLabels;
      for k = 1:length(properties)
         P.other.(properties{k}) = F(obj.other.(properties{k}));
      end
    end
    
    function P = reduce (obj, S)
      % Reduce label data to specified item subset
      % P = obj.reduce(S): Reduce all label data to the item subset M
      %   specified by S. S may be one of the following:
      %     - a boolean vector of length == numItems; is copied to M
      %     - a regular expression name pattern or an index vector selecting
      %       a set of annotations; M is the union of all annotations
      %     - an MSLabelData object; M is the set of all items with non-zero
      %       label values
      %   The resulting new MSDataPartition object is returned.

      P = obj.applyLabelFunction(@(x) x.reduce(S));
    end
    
    function P = removeUnused (obj)
      % Remove unused labels
      % P = obj.removeUnused(): Create a copy of obj with all unused labels
      %   in all label objects removed
      %   The resulting new MSDataPartition object is returned.

      P = obj.applyLabelFunction(@(x) x.removeUnused());
    end
    
    function assert (obj)
      % Check object consistency, abort if inconsistent
      
      % Class labels need to be defined
      assert(~isempty(obj.classes), ...
             'Assertion failed: MSDataPartition.classes must be defined');
      obj.checkLabelProperty(obj.classes, 'classes', []);
      % Check all other properties, need to have equal data size
      % Defined properties
      properties = {'samples', 'individuals', 'sampleSets', 'cvEntities'};
      for k = 1:length(properties)
        obj.checkLabelProperty(obj.(properties{k}), properties{k}, obj.classes.numItems);
      end
      % Fields of struct 'other' must also be properly sized label objects
      properties = fieldnames(obj.other);
      for k = 1:length(properties)
        obj.checkLabelProperty(obj.other.(properties{k}), ...
                               ['other.' properties{k}], obj.classes.numItems);
      end
    end
  end

  methods (Static)
    function L = generatePropertyLabel (D, P, S)
      % Generate label from region property, with optional aggregation
      %
      % L = MSDataPartition.generatePropertyLabel(D, P, S)
      %   D: MS(Maldi)Data object for which label is generated
      %   P: Cell vector of region property values from which unique labels 
      %      and corresponding label values are generated
      %   S: Optional MSLabelData object defining a higher level partition 
      %      over which the generated partition shall be aggregated

      % Find unique property values (= label names) and corresponding map to
      % region indices
      regionIndices = find(~cellfun(@isempty, P));
      [labelNames, ~, regionLabelMap] = unique(P(regionIndices));
      % Generate cell array of regions associated with each label name
      labelRegions = cell(length(labelNames),1);
      for k = 1:length(regionLabelMap)
        j = regionLabelMap(k);
        labelRegions{j} = [labelRegions{j} regionIndices(k)];
      end
      % Generate partition label
      L = D.partition(labelNames, labelRegions);

      % If aggregation label is specified, perform aggregation and inverse
      % aggregation
      if nargin >= 3
        L = L.aggregate(S, 'unique').aggregateInverse(S);
      end
    end
  end
  
  methods (Static, Access = protected)
    function checkLabelProperty (L, propertyName, numItems)
      % Check whether specified property is either empty or an MSLabelData
      % object with one label column. If numItems is specified, it must 
      % match property's numItems.
      if ~isempty(L)
        assert(isa(L, 'MSLabelData'), ...
               'Assertion failed: MSDataPartition.%s must be an MSLabelData object', ...
               propertyName);
        L.assert;
        assert(L.dataLength == 1, ...
               'Assertion failed: MSDataPartition.%s must have a single label column', ...
               propertyName);
        if ~isempty(numItems)
          assert(L.numItems == numItems, ...
                 ['Assertion failed: MSDataPartition.%s.numItems must match ' ...
                  'common number of items (= %d)'], propertyName, numItems);
        end
      end
    end
  end
end

