function mergedT = MSMergeAndTranslateDecStructs(varargin)

  merged = mergeStructs(varargin{:});
  mergedT = MSTranslateDecEquivInputs(merged);

end