function md = MSMaldiUnion (mdList)
% Given a cell array with maldi data objects return their union
% md=MSMaldiUnion({maldiData1,maldiData2,...,maldiDataN})
%   INPUT
%   MdList: cell array of MSMaldiData
%   OUTPUT
%   md: MSMaldiData union

% --------- input validation------------------------------------------
if ~iscell(mdList)||length(mdList)<2
    error('The input must be a cell array of Maldi data of at least two elements');
end
N=length(mdList);
for i=1:N
    if ~isa(mdList{i},'MSMaldiData')
        error('The element %d of the cell array is not an instance of MSMaldiData',i); 
    end
    mdList{i}.assert
end
% ---------- compute common mz interval ---------------------------------
[mzVectorNew, indexFirstCut]=MSCommonMzVector(mdList);
%------- build new data-----------------------------------
disp('Building data union...')
numItemsTotal=0;
allSingle=true;
for i=1:N
    numItemsTotal=numItemsTotal+mdList{i}.numItems;
    allSingle = allSingle & strcmp(class(mdList{i}.data), 'single');
end
if allSingle
  outType='single';
else
  outType='double';
end
data=zeros(numItemsTotal, length(mzVectorNew), outType);
data(1:mdList{1}.numItems,:)=mdList{1}.data(:,indexFirstCut(1):indexFirstCut(2));
updatedPos=mdList{1}.numItems+1;
normalization=mdList{1}.normalization;
for i=2:N
    indexLeft=MSFindMZIndex(mdList{i}.mzVector,mzVectorNew(1));
    indexRight=MSFindMZIndex(mdList{i}.mzVector,mzVectorNew(end));
    % resample only if necessary 
    if (size(mdList{i}.mzVector(indexLeft:indexRight),2)~=size(mzVectorNew,2)||any(mdList{i}.mzVector(indexLeft:indexRight)-mzVectorNew)~=0)
        mdList{i}.resample(mzVectorNew);
    end
    if ~strcmp(mdList{i}.normalization,normalization)
        mdList{i}.setNormalization(normalization);
    end
    data(updatedPos:updatedPos+mdList{i}.numItems-1,:)=mdList{i}.data;
    updatedPos=updatedPos+mdList{i}.numItems;          
end

% --------- Create new maldi data -------------------------------------
md=MSMaldiData(data,mzVectorNew);
disp('Done')
%%-------- build new annotations --------------------------------------
disp('Creating new annotation set...')
annotationSet=annotationsUnion(mdList, numItemsTotal);
md.setAnnotations(annotationSet);
disp('Done')
%--------- build new position grid-------------------------------------
disp('Creating new position grid...')
positionGrid=positionGridUnion(mdList);
if~isempty(positionGrid)
    md.setPositions(positionGrid);
end
disp('Done')
end

function [v, indexFirstCut]=MSCommonMzVector(mdList)
% computes the common mzVector (range intersection with the resolution of
% the first data) from a list of maldi objects.
%  INPUT:
%   mdList: List of maldi objects
%  OUTPUT:
%   v: common mzVector
%   indexFirstCut: two element array with indexes of the original mzVector
%                  of the first element where the min and max common mz
%                  values are found.
%                  
mzLeft=mdList{1}.mzVector(1);
mzRight=mdList{1}.mzVector(end);
N=length(mdList);
for i=2:N
    v=mdList{i}.mzVector;
    if (v(1)>mzRight||v(end)<mzLeft)%It is not possible to find any common mz
        error('The mzVector of the maldi data in position %d has no intersection with one or more of the rest mzVectors',i);
    end
    mzLeft=max(v(1),mzLeft);
    mzRight=min(mzRight, v(end));
end
v=mdList{1}.mzVector;
indexes=v>=mzLeft & v<=mzRight;
v=v(indexes);
indexFirstCut=[find(indexes,1) find(indexes,1,'last')];
end

function annotationSet=annotationsUnion(mdList, numItemsTotal)
N=length(mdList);
if nargin<2
   numItemsTotal=0;
   for i=1:N
       numItemsTotal=numItemsTotal+mdList{i}.numItems;
   end
end

annotationSet=MSAnnotationSet(numItemsTotal);
mdProcessed=0;
for i=1:N
    annot=mdList{i}.annotations;
    if isa(annot,'MSAnnotationSet')
        numAnnot=annot.numAnnotations;
        annotationSet.append({mdList{i}.annotations.annotations(:).name},...
            [false(mdProcessed,numAnnot);...
            mdList{i}.annotations.annotations(:).mask;...
            false(numItemsTotal-mdProcessed-mdList{i}.numItems,numAnnot)]);
    end    
    % add extra annotation indicating the spectra corresponding to
    % independent maldiData
    annotationSet.append({strcat('Union',num2str(i))},[false(mdProcessed,1);...
        true(mdList{i}.numItems,1);false(numItemsTotal-mdProcessed-mdList{i}.numItems,1)]);
    mdProcessed=mdProcessed+mdList{i}.numItems;
end
end

function positionGrid=positionGridUnion(mdList)
% Position grids are placed one under the other
width=0;
% compute width of the indexGrid 
N=length(mdList);
for i=1:N
    pos=mdList{i}.positions;
    if isa(pos,'MSPositionGrid')
        width=max(width,size(pos.indexGrid,2));
    end
end
if width==0 % obtain positions if for any of the given data there was a position representation
    positionGrid=[];
else
    % compute height of the indexGrid 
    height=0;
    for i=1:N
        pos=mdList{i}.positions;
        if isa(pos,'MSPositionGrid')
            height=height+size(pos.indexGrid,1);
        else
            height=height+ceil(mdList{i}.numItems/width);
        end
    end
    mdProcessed=0;
    offset=10; % separation between two maldi data
    indexGrid=zeros(height+(N-1)*offset,width);
    updatedPos=1;
    % ensemble the indexGrids
    for i=1:N
        pos=mdList{i}.positions;        
        if ~isa(pos,'MSPositionGrid')                
            indexG=reshape((1:mdList{i}.numItems), [ceil(mdList{i}.numItems/width),width]);
        else
            indexG=pos.indexGrid;
        end
        indexGrid(updatedPos+offset*(i-1):updatedPos+size(indexG,1)-1+offset*(i-1)...
                                               ,1:size(indexG,2))=indexG+mdProcessed*(indexG>0);

        updatedPos=updatedPos+size(indexG,1);            
        mdProcessed=mdProcessed+mdList{i}.numItems;
        clear indexG;
    end
    positionGrid=MSPositionGrid(indexGrid);
end
end
