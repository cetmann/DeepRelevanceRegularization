classdef MSVal
  % validations functions
  
  methods(Static)
    function bool=isPosInt(x)
      bool=MSVal.isPosNum(x)&&(floor(x)==x);
    end
    
    function bool=isNonNegInt(x)
      bool=MSVal.isPosInt(x)|| (isnumeric(x)&&x==0);
    end
    
    function bool=isPosNum(x)
      bool=isscalar(x)&&isnumeric(x)&&x>0;
    end
  end
  
end
