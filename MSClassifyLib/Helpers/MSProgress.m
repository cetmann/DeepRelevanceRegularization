classdef MSProgress < handle
  % Utility class for progress indicators
  % (See demo method for usage info)
  %
  % Methods:
  %   MSProgress: Constructor
  %   update:     Update progress indicator
  %   close:      Close and end progress indicator
  %   isActive:   Return true if this is the outermost active indicator
  % Static methods:
  %   on:         Globally enable progress indicators
  %   off:        Globally disable progress indicators
  %   clear:      Clear global status
  %   idle:       Return true if no progress indicator is active
  %   demo:       Test and demo code
  %
  
  properties (Access = protected)
    nSteps = 20; % Number of progress steps
    nCurrent;    % Current number of steps
    vMax;        % Maximum progress value
    silent;      % Set to true to hide progress indicator
    ticStart;    % Start time
  end
  
  properties (Dependent)
    visible;     % Visibility of progress indicator
  end
  
  methods
    function obj = MSProgress (name, vMax, silent)
      % Constructor
      %
      % obj = MSProgress(name, vMax): Create progress indicator with given
      %   name and maximum progress value (positive scalar)
      % obj = MSProgress(name, vMax, silent): Specify silent flag (logical
      %   true/false or string 'silent'). Set to true/'silent' to suppress
      %   progress indicator output.
      
      narginchk(2,3);
      
      if ~ischar(name)
        error('Name must be a string');
      end
      if ~(isscalar(vMax) && isnumeric(vMax) && vMax > 0)
        error('Maximum progress value must be a positive scalar');
      end
      if nargin >= 3
        if ischar(silent) && strcmpi(silent, 'silent')
          obj.silent = true;
        elseif isscalar(silent) && islogical(silent)
          obj.silent = silent;
        else
          error('Silent argument must be a logical value or the string ''silent''');
        end
      else
        obj.silent = false;
      end
      
      obj.vMax = vMax;
      obj.reset();
      if obj.idle
        obj.status(obj);
      end
      if obj.visible
        fprintf('%s.', name);
      end
      obj.ticStart = tic;
    end
    
    function update (obj, v)
      % Update progress indicator
      
      if ~(isscalar(v) && isnumeric(v))
        error('Progress value must be a scalar');
      end
      
      n = floor(obj.nSteps*v/obj.vMax);
      if n > obj.nCurrent
        if obj.visible
          fprintf(repmat('.', [1 n-obj.nCurrent]));
        end
        obj.nCurrent = n;
      end
    end
    
    function close (obj)
      % Close and end progress indicator
      
      if obj.visible
        fprintf('done (duration %g sec)\n', toc(obj.ticStart));
      end
      if obj.isActive
        obj.clear();
      end
      obj.reset();
    end
    
    function v = get.visible (obj)
      % Get visibility status
      v = ~obj.silent && obj.isActive;
    end
    
    function a = isActive (obj)
      % Return true if obj is the active MSProgress object
      s = obj.status;
      a = ~isempty(s) && s == obj;
    end
  end
  
  methods (Access = protected)
    function reset (obj)
      % Initialize
      obj.nCurrent = 0;
    end
  end
  
  methods (Static)
    function i = idle
      % Return true if enabled and no progress object active
      i = isempty(MSProgress.status);
    end
    
    function on
      % Globally enable progress indicators
      MSProgress.clear();
    end
    
    function off
      % Globally disable progress indicators
      MSProgress.status(-1);
    end
    
    function clear
      % Clear status
      MSProgress.status([]);
    end
  end
  
  methods (Static, Access = protected)
    function sOut = status (sIn)
      % Access to global status variable
      %
      % sOut = MSProgress.status(): Return current status
      % MSProgress.status(sIn): Set status to sIn
      % sOut = MSProgress.status(sIn): Set status to sIn and return
      %   previous status value
      
      narginchk(0,1);
      
      persistent gStatus;
      if nargout >= 1 || nargin == 0
        sOut = gStatus;
      end
      if nargin >= 1
        gStatus = sIn;
      end
    end
  end
  
  % Demo
  methods (Static)
    function demo
      % MSProgress demo
      disp('Use MSProgress to track progress during a for loop')
      n = 100;
      prg = MSProgress('Running demo loop', n);
      for k = 1:n
        prg.update(k);  % Update on every iteration
        pause(1/n);
      end
      prg.close();  % Operation is done
      
      disp(['Several progress indicators may be nested, but only the ' ...
            'outermost will be visible']);
      prg = MSProgress('Running outer loop', n);
      for k = 1:n
        prg.update(k);
        prgNested = MSProgress('Running inner loop', n);
        for j = 1:n
          prgNested.update(j);
        end
        prgNested.close();
        pause(1/n);
      end
      prg.close();  % Operation is done
      
      disp('Hide a progress indicator by a silent flag in the constructor')
      disp('(No progress indicator visible)');
      prg = MSProgress('Running demo loop (silent) ', n, 'silent');
      for k = 1:n; prg.update(k); pause(1/n); end
      prg.close();
      
      disp('Visibility can also be controlled globally');
      disp('By default, progress indicators are visible');
      prg = MSProgress('Running demo loop', n);
      for k = 1:n; prg.update(k); pause(1/n); end
      prg.close();
      disp('Progress indicators switched off ...')
      MSProgress.off
      prg = MSProgress('Running demo loop', n);
      for k = 1:n; prg.update(k); pause(1/n); end
      prg.close();
      disp('... and switched on again')
      MSProgress.on
      prg = MSProgress('Running demo loop', n);
      for k = 1:n; prg.update(k); pause(1/n); end
      prg.close();
      
    end
  end
end
