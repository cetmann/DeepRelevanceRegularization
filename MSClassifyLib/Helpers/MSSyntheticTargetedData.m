classdef MSSyntheticTargetedData
    % Generates a synthetic data object using targeted proteins. Can be
    % used also to generate general synthetic data (see constructor)
    
    properties
        targetedPeptideIntensities; % cell with the intensities of each targeted protein
        msdata; % generated MSMaldiData
        K; % channel matrix determining the distribution of the 'artificial' proteins
        X; % 'artificial' proteins (row-wise) generating the spectra
        pattern; % MSPatternCell object containing the targeted proteins pattern
        nTargeted; % number of targeted proteins
    end
    methods        
        function obj = MSSyntheticTargetedData(patternCell, numItems, nProteins, noise, nPeptPerProtInt, localizationGridDim) 
            % Generates synthetic targeted data according to the given
            % pattern
            %   INPUT
            %       -patternCell: numeric vector corresponding to an
            %       mzVector (uses no targeted proteins) or an MSPatternCell
            %       object containing the targeted proteins. 
            %       -numItems: number of spectra generated (number of rows
            %       in resulting msdata
            %       nProteins: number of proteins (basis, or X rows) used
            %       to generate data
            %       -noise (optional): a scalar indicating noise level. If
            %       negative the noise is considered multiplicative and
            %       additive for positive numbers.
            %       -nPeptPerProtInt: interval indicating the range for the
            %       number of peptides per protein.
           narginchk(3,6)
           if nargin < 6 || isempty(localizationGridDim)
               localizationGridDim = [5,5];
           elseif nargin < 5 || isempty(nPeptPerProtInt)
               nPeptPerProtInt=[4 6];
           else
               if ~(isnumeric(nPeptPerProtInt)&& size(nPeptPerProtInt,1)==1 &&...
                       size(nPeptPerProtInt,2)==2 && all(nPeptPerProtInt>0) &&...
                       all(nPeptPerProtInt==floor(nPeptPerProtInt)) && ...
                       nPeptPerProtInt(1)<=nPeptPerProtInt(2))
                       
                   error('The 5th input argument must be a valid integer interval')                   
               end
           end
           if nargin < 4 || isempty(noise)
               noise = 0;              
           elseif~(isscalar(noise)&&isnumeric(noise))
               error('the fourth input parameter must be a numerical scalar')
           end
            if isa(patternCell,'MSPatternCell')
                patternCell=reducePatternCell(patternCell,nPeptPerProtInt(1),nPeptPerProtInt(2));
                obj.nTargeted = patternCell.nPatterns;
                mzVector = patternCell.mzVector;
                nExtraProteins = nProteins-obj.nTargeted;
            elseif isnumeric(patternCell)&& size(patternCell,1)==1
                nExtraProteins = nProteins;
                mzVector = patternCell;
                obj.nTargeted = 0;
            else
                error('You must specify either a patternCell object or an mzVector')
            end
            if nExtraProteins < 0
                nExtraProteins = 0;
                nProteins = obj.nTargeted;
                warning(['The number of targeted proteins is larger than the requested proteins'....
                    ' Therefore the minimum number of proteins (as many as targeted) will be generated'])
            end
            % initialize X matrix
            obj.X = zeros(nProteins, length(mzVector));
            
            % initialize K matrix
            
            % Generate chanels randomly
%             obj.K = (rand([numItems, nProteins]));
            
            % Build position grid
            a = sqrt(numItems);
            gridRows = ceil(a);
            gridColumns = ceil(numItems/a);
            positions = reshape([1:numItems zeros(1,gridRows*gridColumns-numItems)],gridRows,gridColumns);
            msPosGrid = MSPositionGrid(positions);
            
            %shuffle seed
            rng('shuffle')
            % build low resolution grid
            zeroProb = 0.6;
            minCoefInterval = [0 5];
            for i=1:nProteins
                obj.K(:,i) = obj.generateChannel([gridRows, gridColumns], localizationGridDim, zeroProb, minCoefInterval, numItems);
            end
            
            % GENERATE TARGETED PROTEINS FROM PATTERN
            if obj.nTargeted                
                obj.targetedPeptideIntensities = cell(1,obj.nTargeted);
                for i=1:obj.nTargeted
                    nPeptides = size(patternCell.patterns{i},1);
                    obj.targetedPeptideIntensities{i} = (1 + rand([1,nPeptides]))/2;
                    obj.X(i,:) = obj.X(i,:) + obj.targetedPeptideIntensities{i}*patternCell.patterns{i};
                    obj.X(i,:) = obj.X(i,:)/norm(obj.X(i,:));
                end
            end

            % GENERATE NON-TARGETED RANDOM PROTEINS
            
            % obtain number of peptides for the non-targeted proteins
            k = randi(nPeptPerProtInt,1,nExtraProteins);

            for j=1:nExtraProteins
                % obtain the mz peaks for the current synthetic protein
                peaks = datasample(ceil(mzVector(1)):floor(mzVector(end-20)),k(j),'Replace',false);
                for i=1:k(j)
                    % obtain the isotopic pattern corresponding to the peak
                    [~,isoPattern] = isotopicPattern(peaks(i)*1.0005,mzVector,0.2,0);
                    % obtain the intensity of this peptide on the synthetic
                    % data
                    intensity = (rand);
                    obj.X(obj.nTargeted+j,:)=...
                        obj.X(obj.nTargeted+j,:)+...
                        isoPattern*intensity;
                end
                obj.X(obj.nTargeted + j,:)=obj.X(obj.nTargeted + j,:)/(norm(obj.X(obj.nTargeted + j,:)));
            end
            % initialize data
            data = obj.K*obj.X;
            
            if noise > 0
                % Add additive noise to data
                data = data + noise*randn(size(data));
                data(data<0)=0;
            elseif noise<0            
                % Add multiplicative noise to data            
                data = data.*(abs(1 + noise*randn(size(data))));
            end
            
            % create maldi synthetic data              
            obj.msdata = MSMaldiData(data, mzVector);
            obj.msdata.setPositions(msPosGrid);
            obj.pattern = patternCell;
        end
    end
    methods (Static)
        function v = generateChannel(positionsDim, localizationGridDim, zeroProb, minCoefInterval, nElements)
            v = zeros(positionsDim);
            localizationGrid = rand(localizationGridDim);
            localizationGrid = localizationGrid > zeroProb;
            localizationGrid = localizationGrid.*...
                  (minCoefInterval(1)+(minCoefInterval(2)-minCoefInterval(1))*rand(localizationGridDim));
            nR = positionsDim(1);
            nC = positionsDim(2);
            hStep = ceil(nC/localizationGridDim(2));
            vStep = ceil(nR/localizationGridDim(1));
            for i=1:localizationGridDim(1)
                initV = (i-1)*vStep+1;
                endV = min(i*vStep,nR);
                for j=1:localizationGridDim(2)
                    if localizationGrid(i,j)
                        initH = (j-1)*hStep+1;
                        endH = min(j*hStep, nC);
                        v(initV:endV,initH:endH) = localizationGrid(i,j)+(localizationGrid(i,j)+1)*...
                                                   rand(endV-initV+1,endH-initH+1);
                    end
                end
            end
           v = v(1:nElements);
        end
    end
end