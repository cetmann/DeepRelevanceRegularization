function merged = mergeStructs(varargin)
p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});
merged = p.Unmatched;
end
