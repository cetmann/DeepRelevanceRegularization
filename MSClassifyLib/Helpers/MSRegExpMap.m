function newStr = MSRegExpMap(str, expression, replace, varargin)
  % Map string values using regular expression patterns
  %
  % newStr = MSRegExpMap(str, expression, replace, ...)
  %   str: Cell array of input strings
  %   expression: Regular expression search pattern
  %   replace: Regular expression replacement pattern
  %   ...: Additional arguments passed to regexp(rep) (optional)
  %   newStr: Cell array of replacement strings, empty if no match
  
  expression = ['^.*?' expression '.*?$'];
  newStr = cell(size(str));
  inMask = ~cellfun(@isempty, str);
  matchMask = ~cellfun(@isempty, regexp(str(inMask), expression, 'once', varargin{:}));
  outMask = inMask;
  outMask(outMask) = matchMask;
  newStr(outMask) = regexprep(str(outMask), expression, replace, varargin{:});
end
