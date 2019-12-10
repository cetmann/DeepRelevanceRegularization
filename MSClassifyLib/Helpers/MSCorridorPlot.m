function varargout = MSCorridorPlot (X, Y, L, U, varargin)
  % Draw a corridor plot consisting of a center line and an error corridor
  % h = MSCorridorPlot (X, Y, L, U, ...)
  %   Parameters are as in errorplot
  %
  % TODO: More detailed description

  % Dummy plot center line to determine its color
  % Retain current axes hold status
  nextPlotSave = get(gca, 'NextPlot');
  % If next plot will add to existing plots, retain current color order index
  if strcmp(nextPlotSave, 'add')
    colorOrderIndexSave = get(gca, 'ColorOrderIndex');
  else
    % Otherwise assume that color order index will be reset to 1
    colorOrderIndexSave = 1;
  end
  % Plot center line and obtain line color
  hLine = plot(X, Y, varargin{:});
  lineColor = get(hLine, 'Color');
  % Remove line and restore previous color order index
  delete(hLine);
  set(gca, 'ColorOrderIndex', colorOrderIndexSave);

  % Retain current axes hold status and temporarily set to hold on
  nextPlotSave = get(gca, 'NextPlot');
  hold on

  % Avoid plotting the corridor in case the bounds are all zero, as the
  % Matlab fill() function does not seem to properly handle this case.
  if ~(any(L(:)) || any(U(:)))
    hCorridor = 0;
  else
    % Plot corridor as filled, transparent polygon without outline
    hCorridor = fill([X; flipud(X)], [Y-L; flipud(Y+U)], lineColor, ...
                     'LineStyle', 'none');
    set(hCorridor, 'FaceAlpha', 0.1);
  end
  
  % Restore color index
  set(gca, 'ColorOrderIndex', colorOrderIndexSave);
  % Plot center line on top of corridor polygon
  hLine = plot(X, Y, varargin{:});

  % Restore hold status
  set(gca, 'NextPlot', nextPlotSave);

  % If output argument is specified, return vector of object handles
  if nargout >= 1
    varargout{1} = [hCorridor; hLine];
  end

end

