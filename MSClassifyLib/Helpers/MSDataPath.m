function path = MSDataPath (varargin)
  % Set or get MSClassifyLib data search path
  % The data search path is a list of directories stored as a cell string
  % array. It is used by MSResolvePath to locate a data file or directory,
  % with directories listed first having precedence over those listed
  % subsequently.
  % path = MSDataPath(): Return current data search path
  % path = MSDataPath(p1, p2, ...): Set data search path to the list p1,
  %   p2, etc. Each argument pk may be a string or a cell string array. The
  %   resulting new search path is returned.
  % path = MSDataPath([]): Clear data search path
  
  path = {};
  if nargin == 0
    % No argument specified, retrieve current data search path
    if ispref('MSClassifyLib', 'dataPath')
      path = getpref('MSClassifyLib', 'dataPath');
    end
    
  else
    % Arguments specify new data search path
    for k = 1:length(varargin)
      % Argument may be either a string or a cell string array
      if ischar(varargin{k})
        path = [path; varargin(k)];
      elseif iscellstr(varargin{k})
        path = [path; varargin{k}(:)];
      elseif ~isempty(varargin{k})
        error('Arguments must be either strings or cell string arrays');
      end
    end
    
    if isempty(path)
      % Empty path specified, clear search path
      rmpref('MSClassifyLib', 'dataPath');
    else
      % Remove duplicate directories
      path = unique(path, 'stable');
      % Store new path to preferences
      setpref('MSClassifyLib', 'dataPath', path);
    end
  end
end

