function P = MSAdaptiveResampleParams (varargin)
  % Parse adaptive resampling parameters and return parameter struct
  % P = MSAdaptiveResampleParams(name, value, ...): Specify parameters for
  %   adaptive resampling as name value pairs. Supported parameters are:
  %   - mzRange: Explicit mass range [mzMin, mzMax] (default = [])
  %   - mzType: 'exact', 'nominal', 'center' (default = 'exact')
  %   - width: Interval width in Da (default = 0.4)
  %   - subsampling: Subsampling factor (positive integer) specifying the
  %     number of subdivisions of 1 Da intervals for the resulting m/z
  %     vector (default = 1)
  %   - massScaleDelta: Relative mass defect shift per mass unit
  %     (default = 4.95e-4, specific to peptides)
  %   - unwrap: Specify whether phase unwrapping is applied (default: false)

  % Validation functions
  isValidRange = @(x) isempty(x) || (isnumeric(x) && isvector(x) && length(x) == 2);
  isValidType = @(x) any(strcmpi(x, {'exact', 'nominal', 'center'}));
  isScalar = @(x) isnumeric(x) && isscalar(x);
  isPosScalar = @(x) isScalar(x) && x > 0;
  isPosInteger = @(x) isPosScalar(x) && mod(x,1) == 0;
  isLogicalScalar = @(x) isscalar(x) && islogical(x);
  % Setup input parser
  params = inputParser;
  params.addParameter('mzRange', [], isValidRange);
  params.addParameter('mzType', 'exact', isValidType);
  params.addParameter('width', 0.4, isPosScalar);
  params.addParameter('subsampling', 1, isPosInteger);
  params.addParameter('massScaleDelta', nan, isScalar);
  params.addParameter('unwrap', false, isLogicalScalar);
  % Parse arguments
  params.parse(varargin{:});
  P = params.Results;
end

