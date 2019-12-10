classdef MSPatternCell
    % Class representing a collection of targeted proteins
    properties(SetAccess=immutable)
        patterns; % cell array. Each element corresponds to a protein and consists of 
                  % a matrix which contains a peptide (including isotopic pattern) row-wise.
        mzVector; % mz-vector
        proteinNames; % optional names of the targeted proteins
    end
    properties(Dependent)
        nPatterns;
    end
    methods
        function obj = MSPatternCell(patterns, mzVector, proteinNames)
            % Constructor
            %   INPUT:
            %       -patterns: cell array. Each element corresponds to a protein and consists of 
            %       a matrix which contains a peptide (including isotopic pattern) row-wise.
            %       -mzVector: mz-vector
            %       -proteinNames: optional names of the targeted proteins
            narginchk(2,3);
            if nargin==3
                if ~iscellstr(proteinNames)||(length(proteinNames)~=length(patterns))
                    error(['proteinNames input must be a cell of strings with '... 
                        'the same length than the input patterns'])
                end
                obj.proteinNames=proteinNames;
            else
                obj.proteinNames = [];
            end
            obj.mzVector = mzVector;
            if ischar(patterns)
                % Load Pattern given by excel table
                try
                    obj.patterns = obj.generatePatternsFromExcel(patterns, mzVector);
                catch e
                    error(['There was a problem generating pattern from table. Error'...
                          'thrown: ' e.message])
                end
            else
                obj.validatePatternCell(patterns);
                obj.patterns = patterns;
            end
        end
        function n = get.nPatterns(obj)
            n = numel(obj.patterns);
        end
    end
    methods (Static)
        function validatePatternCell(patterns)
            % To be implemented
        end
        function patterns = generatePatternsFromExcel(path, mzVector)
            % To be implemented
            patterns = [];
        end
    end
end