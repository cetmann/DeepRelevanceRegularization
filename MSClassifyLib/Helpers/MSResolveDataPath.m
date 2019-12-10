function path = MSResolveDataPath (path, searchPath)
  % Resolve path by searching over current data search path
  %
  % MSResolveDataPath(path): 
  % - If path is relative, it is tested whether it specifies an existing 
  %   file or directory under any of the directories listed in the current 
  %   data search path. If a match is found, this is returned, otherwise
  %   path is returned unmodified. 
  % - If path is absolute, it is returned unmodified.
  % - The data search path may be set or queried using MSDataPath().
  %
  % MSResolveDataPath(path, searchPath): Explicitly specify search path as
  %   a cell string array, ignore search path set by MSDataPath().
  % Update: - Also allowing to specify a relative path in order to return
  % the full path
  % - Returning a warning if the file or folder couldn't be found
  % (otherwise after specifying a non-existing directory errors can occurr
  % much later, possibly difficult to track)

  if ~isRelative(path)
    flagFound = exist(path, 'file') > 0;
  else
    if exist(path, 'file')
      flagFound = true;
    else
      if nargin < 2
        searchPath = MSDataPath;
      end
      flagFound = false;
      for k = 1:length(searchPath)
        % convert all slashes to usual slashes of working system  
        if strcmp('/', filesep)  
            searchPath{k} = strrep(searchPath{k}, '\', filesep);
            path = strrep(path, '\', filesep);
        else
            searchPath{k} = strrep(searchPath{k}, '/', filesep);
            path = strrep(path, '/', filesep);
        end
        p = fullfile(searchPath{k}, path);
        if exist(p, 'file')
          flagFound = true;
          path = p;
          break;
        end
      end
    end
  end
  if ~flagFound
    warning('File / folder not found');
  end
end

function R = isRelative (path)
  % Return true if path is relative on the current platform
  
  if ispc
    % On Windows, an absolute path starts with a slash (forward or 
    % backward), optionally preceded by a drive letter and colon.
    R = isempty(regexp(path, '^([a-z]:|)[\\/]', 'ignorecase', 'once'));
  elseif isunix
    % On Unix, an absolute path starts with '/' or '~/'
    R = isempty(regexp(path, '^(~|)/', 'once'));
  else
    error('isRelative() not implemented for this platform')
  end
end
