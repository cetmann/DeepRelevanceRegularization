classdef slDump < handle
    % Read access to Scilslab(R) .sl files.
    %
    % This class dumps selected parts of sl files:
    % * region information (Names, Properties, Spots)
    % * spectra-matrix
    %
    % Properties
    %    RegionNames        : nx1 cell : full name of the regions.     (data)
    %    RegionSpots        : nx1 cell : spot indices for each region  (data)
    %    RegionProperties   : nxm cell : region properties             (data)
    %    RegionPropertyNames: 1xm cell : name of region properties     (data)
        %    
    %    raster             : struct   : parameters of rasterization    (dev)
    %    iscompressed       : 0|1      : switch for LZ4 compression    (info)
    %    version            : int      : version number of sl scheme   (info)
    %
    % Methods:
    %    getSpectraMatrix   : returns the specra matrix
    %    getRegionTable     : returns table with collected region information
    %    selftest           : test for knwon bugs / problems (should be
    %                         executed once before processing new datasets)
    %                         return value must be 1 (no errors)
    %
    % slDump uses the handle semantic, i.e. when assigning an object
    % of this class to a variable, only a reference to the original object is
    % copied. 
    
    properties(SetAccess=protected)
        filename;
        
        RegionNames={};
        RegionSpots={};
        RegionProperties={}; 
        RegionPropertyNames={};
        
        iscompressed=0;
        
        version;
               
    end
    
    properties
        verbose=0;
        debug=0;
        raster;
    end
    
    properties(SetAccess=protected,GetAccess=protected)
        urlRegionDataSet={};
        rootreg=1;
        PropMapInt; % region property hexid -> Integer (column index)
                    % for regionproperties (deprecated)  
        
        AttributeRegions; %internal data structure with names, types and urls
    end
    
    methods
        
        function this=slDump(filename)
            slDump.checkplugin();
            if ~exist(filename,'file'),
                error('File not found (%s)', filename);
            end
            this.filename=filename;
            this.version=h5read(filename,'/Version');
            if this.version>14,
                cflag=h5readatt(filename,'/','Compressed');
                if(strcmp(cflag,'true')),
                    this.iscompressed=1;
                    if verLessThan('matlab','8.5'),
                        warning('matlab >= 8.5 (2015a) is required.');
                    end
                end
            end
            if this.iscompressed,  c='c';else c='';  end
            fprintf('loading %s (v%d%s)\n', this.filename, this.version, c);
            
            %fprintf('scanning region property table ...')
            %this.initRegionProperties;
            %fprintf('ok (%d region properties found)\n', numel(this.RegionPropertyNames));
            
            fprintf('scanning regions ...');
            this.initregions();
            fprintf('ok (%d regions found)\n', numel(this.urlRegionDataSet));
            
            
            fprintf('scanning regions spots ...');
            this.RegionNames=cell(size(this.urlRegionDataSet));
            this.RegionSpots=cell(size(this.urlRegionDataSet));
            for i=1:length(this.urlRegionDataSet)
                this.RegionNames{i}=this.recH5queryName(this.urlRegionDataSet{i});
                this.RegionSpots{i}=h5read(this.filename, [this.urlRegionDataSet{i},'/SpotList'])+1;
            end
            this.RegionNames=this.RegionNames';
            this.RegionSpots=this.RegionSpots';
            %this.RegionMetaData=this.RegionMetaData';
            fprintf(' ok\n');
            
            this.raster.rasteroffset=[];
            this.raster.rasterwhitespacescaling=[];
            this.raster.rasterreg=[];
            this.raster.rastersize=[];
            
            
            
            fprintf('scanning Attributes ...');
            H=this.slInfo('/AttributeRegions', [],31);
            if ~isempty(H)
                urls=this.reccollectUrls(H.Groups, 'AttributeRegion');
                fprintf(' ok (%d found)\n', numel(urls));
            else
                urls=[];
                fprintf(' skipped\n');
            end
                                  
            this.AttributeRegions.urls=urls;
           
            if ~isempty(urls)
                fprintf('Read RegionProperties ...');
                for i=1:numel(urls)
                    this.AttributeRegions.name{i}=this.slreadatt(urls{i}, 'name');
                    this.AttributeRegions.ValueType(i)=this.slread([urls{i},'/ValueType']);
                end
                this.initAttributeTable();
                fprintf(' ok \n');
            else
                this.AttributeRegions.name={};
                this.AttributeRegions.ValueType=[];
            end
        end
        
        function A=getSpectra(this, I, J_range)
            % Load spectral data for spot indices I and m/z index range 
            % J_range (optional, defaults to full m/z range).
            
            if nargin < 3
                % m/z range is optional
                J_range=[];
            end
            
            % Convert index set I to list of intervals
            iList=slDump.indexset2intervals(I);
            % Check spectra index limit
            nspec=this.getNumberOfSpectra();
            if max(iList(:,2)) > nspec
                error('out of bounds; %d spectra are available.', nspec);
            end
            
            % Check m/z range
            nmz=length(this.getMZvec);
            isrange=@(x) isnumeric(x) && isvector(x) && length(x)==2 && ...
                         x(1)>=1 && x(2)>=x(1) && all(mod(x,1)==0);
            if isempty(J_range)
                J_range=[1 nmz];
            elseif ~isrange(J_range)
                error('invalid m/z index range argument');
            elseif J_range(2)>nmz
                error('out of bounds; %d m/z values are available.', nmz);
            end
            numJ = diff(J_range)+1;
            
            if(this.debug>0),
                h5disp(this.filename, '/SpectraGroups/InitialMeasurement/spectra');
            end
            
            % Load data
            A=nan(numJ,length(I),'single');
            row=1;
            for k=1:size(iList,1)
                numI=iList(k,2)-iList(k,1)+1;
                A(:,row:row+numI-1)=...
                  h5read(this.filename, '/SpectraGroups/InitialMeasurement/spectra', ...
                         [J_range(1) iList(k,1)],[numJ numI]);
                row=row+numI;
            end
        end
        
        function S=getSpectraMatrix(this)
            S=h5read(this.filename, '/SpectraGroups/InitialMeasurement/spectra');
        end
              
        function A=getMZdata(this, I)
            [istart,iend]=slDump.indexset2interval(I);
            
            num=iend-istart+1;
            nspec=this.getNumberOfSpectra();
            nmz=length(this.getMZvec);
            
            if iend > nmz
                error('out of bounds; %d mz values are available.',nmz);
            end            
            if(this.debug>0),
                h5disp(this.filename, '/SpectraGroups/InitialMeasurement/images');
            end            
            A=h5read(this.filename, '/SpectraGroups/InitialMeasurement/images',[1 istart],[nspec num]);
        end
        
        function S=getMZdataMatrix(this)
            S=h5read(this.filename, '/SpectraGroups/InitialMeasurement/images');
        end
              
        
        function mz=getMZvec(this)
            mz=h5read(this.filename,'/SamplePositions/GlobalMassAxis/SamplePositions');
            mz=double(mz);
        end
        
        function M=getMZimagesI(this,I)
            A=this.getMZdata(I);
            [istart,iend]=slDump.indexset2interval(I);
            IDX=this.getIndexMatrix();
            
            M=nan(size(IDX,1)*size(IDX,2), iend-istart+1, 'like', A);
            
            M(IDX~=0,:)=A;
            M=reshape(M,size(IDX,1), size(IDX,2),[]);
        end

        function i=getRasterIndices(this)
            i=h5read(this.filename,'/MeasurementRasterIndices');
        end
        
        function T=getRasterPixel2WorldTrans(this, hexid)
            % hexid for old h5 files needed
            if nargin < 2,
                hexid='0';
            end
            
            I=this.getRasterIndices();
            uniqi=unique(I);
            T=zeros(4,4,length(uniqi));
            for i=1:length(uniqi),
                T(:,:,i)=h5read(this.filename, sprintf('/Registrations/%s/Rasters/%d/PixelToWorld', hexid, uniqi(i)))';
            end
        end
        
        function C=getWorldCoordinates(this, hexid)
            %hexid for old h5 files needed
            if nargin < 2,
                hexid='0';
            end
            C=h5read(this.filename, sprintf('/Registrations/%s/Coordinates', hexid));
            C(:,4)=1;
            C=C';
        end
        
        function C=getCoordinates(this, hexid)
            warning('Please use "getWorldCoordinates" instead of "getCoordinates" in future.');
            if nargin < 2,
                C=this.getWorldCoordinates();
            else
                C=this.getWorldCoordinates(hexid);
            end
        end
        
        
        
        function T=getRegionTable(this) %deprecated
            if isempty(this.RegionProperties),
                T=table();
            else
                %find all empty cols
                I=all(cellfun(@isempty,this.RegionProperties));
                %ensure valid and unique identifiers for table
                [validStrings,      modV]=matlab.lang.makeValidName(this.RegionPropertyNames);
                [validUniqueStrings,modU]=matlab.lang.makeUniqueStrings(validStrings,{}, namelengthmax);
                T=cell2table(this.RegionProperties, 'VariableNames',validUniqueStrings);
                
                for i=find(modV|modU),
                    warning('Redefinition of Property name: %s => %s', this.RegionPropertyNames{i}, validUniqueStrings{i});
                end
                
                if any(I),
                    %remove empty cols
                    T(:,I)=[];
                    warning('Ignored %d empty columns.', sum(I));
                end
            end
            %concat Region Names
            T=[cell2table(this.RegionNames, 'VariableNames',{'RegionName'}) ,T];
            %concat nSpots
            T=[T,table(cellfun(@numel, this.RegionSpots), 'VariableNames',{'nSpots'})];
            %row names 1:n
            num=cell(size(this.RegionNames,1),1);
            for i=1:length(num), num{i}=num2str(i); end
            T.Properties.RowNames=num;
        end
        
        function enlargewhitespace(this)
            
            this.raster.rasteroffset=[];
            this.raster.rasterwhitespacescaling = this.raster.rasterwhitespacescaling / 2;
            this.raster.rasterreg=[];
            this.raster.rastersize=[];
            
        end
        
        
        
        function R = getRasterCoordinates(this, r)
            % returns the Raster
            if nargin<2,
                r=1;
            end
            if r~=1,
                error('rasterindex neq 1: to be implemented')
            end
            
            %check state object (global properties must be calculated from
            %root region)
            if(r~=this.raster.rasterreg)
                this.raster.rasterwhitespacescaling=[];
                this.raster.rasteroffset=[];
            end
            
            
            if(~this.rootreg && ...
                    (isempty(this.raster.rasterwhitespacescaling) || isempty(this.raster.rasteroffset)))
                error('Object is in wrong state, getRasterCoordinates() must initially be called for the root region.');
            end
            
            P=this.getWorldCoordinates();
            core=this.getRasterIndices();
            if(size(core,2) > 1)
                warning('%d Rasters available (???) => using no 1', size(core,2));
                core=core(:,1);
            end
            core=core+1;
            uncore=unique(core);
            
            T=this.getRasterPixel2WorldTrans;
            
            %T invertieren
            Tinv=zeros(size(T));
            for i=1:size(T,3),
                Tinv(:,:,i) = T(:,:,i)^(-1);
            end
            
            R=zeros(size(P));
            for i=1:size(R,2),
                R(:,i) = Tinv(:,:,uncore==core(i))*P(:,i);
            end
            
            
            R=round(R);
            
            %offset in world coordinates
            Off=zeros(4, size(P,2));
            for i=1:size(P,2),
                Off(:,i)=T(:,:,uncore==core(i))*[0 0 0 1]';
            end
            
            %calc whitespacescaling factor (once, initially)
            %min bounding box for all cores
            k=this.raster.rasterwhitespacescaling;
            if(isempty(k))
                k=[inf ,inf, inf];
                for i=1:length(uncore),
                    idx=find(core==uncore(i));
                    k=min(k,[
                        (max(P(1,idx)) - min(P(1,idx)))/(max(R(1,idx)) - min(R(1,idx))) ,...
                        (max(P(2,idx)) - min(P(2,idx)))/(max(R(2,idx)) - min(R(2,idx))) ,...
                        (max(P(3,idx)) - min(P(3,idx)))/(max(R(3,idx)) - min(R(3,idx)))
                        ]);
                    
                end
                k=k/sqrt(2); %workaround
                k(1:2)=[min(k) min(k)];
                this.raster.rasterwhitespacescaling=k;
                this.raster.rasterreg=r;
            end
            
            
            %apply scaled offset to raster coords (idea is to remove white space
            %around cores)
            for i=1:3,
                if(k(i) ~= inf),
                    R(i,:)=R(i,:)+Off(i,:)/k(i);
                end
            end
            
            %positive integers ard needed for raster coords
            R=round(R);
            if(isempty(this.raster.rasteroffset)),
                this.raster.rasteroffset=-min(R,[],2)+1;
            end
            for i=1:3,
                R(i,:)=R(i,:)+this.raster.rasteroffset(i);
            end
            
            if isempty(this.raster.rastersize),
                this.raster.rastersize=[max(R(2,:)), max(R(1,:))];
            end
            
        end
        
        function P=getIndexMatrix(this,r)
            if nargin < 2,
                r=1;
            end
            
            
            R=this.getRasterCoordinates(r);
            
            P=zeros(this.raster.rastersize);
            P(mod(R(1,:)-1, size(P,2))*size(P,1) +  R(2,:))=1:size(R,2);
            
        end
        
        
        function R = getRasterCoordinatesRaw(this, r)
            % returns each raster
            if nargin<2,
                r=1;
            end
            if r~=1,
                error('rasterindex neq 1: to be implemented')
            end
            
            P=this.getWorldCoordinates();
            core=this.getRasterIndices();
            if(size(core,2) > 1)
                warning('%d Rasters available (???) => using no 1', size(core,2));
                core=core(:,1);
            end
            core=core+1;
            uncore=unique(core);
            
            T=this.getRasterPixel2WorldTrans;
            
            %T invertieren
            Tinv=zeros(size(T));
            for i=1:size(T,3),
                Tinv(:,:,i) = T(:,:,i)^(-1);
            end
            
            R=zeros(size(P));
            for i=1:size(R,2),
                R(:,i) = Tinv(:,:,uncore==core(i))*P(:,i);
            end
            
            %R=round(R);
            
            
            
        end
        
        
        function P=getIndexMatrixRaw(this,r)
            if nargin < 2,
                r=1;
            end
            
            
            R=this.getRasterCoordinatesRaw(r);
            
            R=round(R);
            
            %calculate overall rastersize
            rastersize=[max(R(2,:)), max(R(1,:))];
            
            core=this.getRasterIndices;
            core=core(:,1);
            core=double(core);
            core=core';
            
            P=zeros(rastersize(1), rastersize(2), length(unique(core)));
            n=size(P);
            
            P(mod(R(1,:), n(2))*n(1) +  R(2,:) +1 +core*(n(1)*n(2)))=1:size(R,2);
            
        end
        
        
        
        function i=getNumberOfSpectra(this)
            i=str2double((cell2mat(h5readatt(this.filename,'/','nSpots'))));
        end
        
        function L=getLabelsH5(this, hexid)
            warning('this is a dev feature');
            L=h5read(this.filename, sprintf('/%s/Labels', hexid));
        end
        
        function err=selftest_raster(this)
            err=0;
            
            % Test the raster
            
            I=this.getIndexMatrix();
            idx=find(I);
            if(numel(unique(I(idx)))~=this.getNumberOfSpectra()),
                
                err=1;
                fprintf('unique positions in raster matrix: %d\n', numel(unique(I(idx))));
                fprintf('Number of Spectra                : %d\n', this.getNumberOfSpectra());
                idx2=find(  ~ismember(1:this.getNumberOfSpectra() , I(idx) ));
                fprintf('No raster point for spectrum     : %d\n', idx2);
                warning('Problem with index matrix');
            end
            
            
        end
        
        
        
        function I=lsAttributes(this, regexstr, casesensitive)
            if nargin<3
                casesensitive=0;
            end
            if nargin<2
                regexstr='.';
            end
            
            if casesensitive
                I=cellfun(@(x) ~isempty(regexp(x,regexstr, 'once')), this.AttributeRegions.name);
            else
                I=cellfun(@(x)~isempty(regexpi(x,regexstr, 'once')), this.AttributeRegions.name);
            end
            
            if this.verbose>=0
                fprintf('%s\n', this.AttributeRegions.name{I})
            end
        end
        
        
        function err=selftest_allregions(this)
            allspots=[];
            for i=1:length(this.RegionSpots),
                allspots=[allspots; this.RegionSpots{i}];
            end
            
            if length( unique(allspots) ) == this.getNumberOfSpectra(),
                err=0;
            else
                err=1;
                fprintf('number of unique spots over all regions = %d\n', length( unique(allspots) ));
                fprintf('total number of spectra= %d\n', this.getNumberOfSpectra() );
                warning('number of spots in regions correct');
            end
        end
        
        
        function ok=selftest(this)
            
            
            err=[];
            
            err(1)=this.selftest_raster();
            err(2)=this.selftest_allregions();
            
            fprintf('%d errors occured\n', sum(err));
            
            
            ok=sum(err)==0;
        end
        
        
        function showinfo (obj)
            % Display information on this object
            
            slInfo = struct;
            [~,name,ext] = fileparts(obj.filename);
            slInfo.file = [name ext];
            slInfo.fullpath = obj.filename;
            slInfo.numberOfSpectra = obj.getNumberOfSpectra();
            mzv = obj.getMZvec()';
            slInfo.dataLength = length(mzv);
            slInfo.mzRange = mzv([1 end]);
            slInfo.mzBinWidthMean = mzv([1 end])*[-1 1]'/slInfo.dataLength;
            slInfo.mzBinWidthRange = mzv([1 2 end-1 end])*[-1 1 0 0; 0 0 -1 1]';
            slInfo.numberOfRegions = length(obj.RegionNames);
            disp(slInfo);
        end
        
    end
    
    methods (Access = protected)
        
        %extract Region Property Names from json string
        function json=initRegionProperties(this)
            try
                json=cell2mat(h5read(this.filename,'/RegionProperties'));
            catch %ME
                warning('missing RegionProperties, file too old?');
                json='';
                %rethrow(ME)
            end
            
            %rudimentary json parser
            [a,b]=regexp(json,'"propertyIds":{[^}]*');
            
            subjson=json(a + length('"propertyIds":{') :  b);
            
            if isempty(subjson),
                this.RegionPropertyNames={};
            else
                splitStr = regexp(subjson,',','split');
                KEY=cell(0);
                VAL=cell(0);
                for i=1:length(splitStr),
                    kv=regexp(splitStr{i},':','split');
                    key=kv{1};
                    val=kv{2};
                    KEY{i}=key(2:end-1);
                    VAL{i}=val(2:end-1);
                end
                
                %set region property names
                this.RegionPropertyNames=VAL;
                
                if( length(unique(VAL)) ~= length(VAL) ),
                    warning('Reading the Region Properties is not possible. Propertiy Names must be unique.');
                end
                
                % map prop hex key -> index
                this.PropMapInt   = containers.Map(KEY, 1:length(KEY));
            end
        end
        
        
        function name=recH5queryName(this, url)
            
            [urlpath,~] = fileparts(url);
            name=cell2mat(h5readatt(this.filename, url, 'name'));
            
            if strcmp(urlpath,'/Regions') || strcmp(urlpath,'/'),
                name=strcat(urlpath(2:end),'/',name);
            else
                name=strcat(this.recH5queryName(urlpath),'/',name);
            end
        end
        
        
        
        
        function counter=regionrec(this, Groups, counter)
            if nargin < 3,
                counter=0;
            end
            
            for i=1:length(Groups)
                %fprintf('%s\n', Groups(i).Name);
                
                %define empty structs to store collected attributes
                %metadata=struct; %empty user metadatax (DEACTIVATED)
                Prop=struct;
                Prop.key=cell(0);
                Prop.val=cell(0);
                
                isregion=false;
                
                for j=1:length(Groups(i).Attributes)
                    key=Groups(i).Attributes(j).Name;
                    val=Groups(i).Attributes(j).Value{1};
                    %fprintf('%15s  =>  %s\n', key , val  );
                    
                    if( strcmp(key, 'type'))
                        if(strcmp(val,'Region')),
                            counter=counter+1;
                            isregion=true;
                            this.urlRegionDataSet{counter}=Groups(i).Name;
                        end
                        
                    elseif( regexp(key, 'User/.*') == 1)
                        %**  User MetaData (Ignored)
                        %fprintf('%15s  =>  %s\n', key , val  );
                        
                        %                       metadata.(key(6:end)) = val;
                        %elseif ( strcmp(key, 'name')),
                        %    metadata.(key) = val;
                        
                        % Check RegionProperty Element
                        if( regexp(key, 'User/property/.*') == 1) % property (table element)
                            Prop.key{length(Prop.key)+1} = key(length('User/property/')+1:end);
                            Prop.val{length(Prop.val)+1} = val;
                        end
                    end
                end
                
                %store collected attributes (RegionProperties & metadata (deactivated))
                if (isregion),
                    %%%this.RegionMetaData{counter} = metadata;
                    for j=1:length(Prop.key),
                        
                        if ~this.PropMapInt.isKey(Prop.key{j})
                            warning('Unnamed region property %s detected, adding', Prop.key{j});
                            % Add new property, using key as name
                            propNameIndex=length(this.RegionPropertyNames)+1;
                            this.RegionPropertyNames{propNameIndex}=Prop.key{j};
                            this.RegionProperties(:, propNameIndex)=cell(size(this.RegionProperties, 1), 1);
                            this.PropMapInt(Prop.key{j})=propNameIndex;
                        end
                        this.RegionProperties{counter, this.PropMapInt(Prop.key{j})}=Prop.val{j};
                    end
                end
                
                %recursion
                if isfield(Groups(i), 'Groups'),
                    counter=this.regionrec(Groups(i).Groups, counter);
                end
                
                
            end
            
            
        end
        function initregions(this)
            %get struct with h5 content for /regions directory
            H=h5info(this.filename,'/Regions');
            %allocate memory for cell arrays
            n=length(H.Groups); %upper bound for number of regions
            this.urlRegionDataSet=cell(n, 1);
            this.urlRegionDataSet=cell(0);
            if ~isempty(this.PropMapInt),
                this.RegionProperties=cell(n, length(this.PropMapInt.values));
                this.RegionProperties=cell(0, length(this.PropMapInt.values));
            end
            %start rec travel
            this.regionrec(H.Groups);
            %fill RegionProperties with empty cells ( urlRegionDataSet has correct length )
            if ~isempty(this.RegionProperties)
                if size(this.RegionProperties, 1) < length(this.urlRegionDataSet),
                    this.RegionProperties{length(this.urlRegionDataSet),length(this.PropMapInt.values)}=[];
                end
            end
            
            
            
        end
        
    %fill the attribute table (region property table)
        function initAttributeTable(this) 
            this.RegionPropertyNames=this.AttributeRegions.name; %take all attributes
            this.RegionProperties=cell(length(this.urlRegionDataSet), length(this.RegionPropertyNames));
            
            for i=1:size(this.RegionProperties,2)
                if this.verbose>0
                    fprintf('%d: %s\n', i, this.AttributeRegions.name{i});
                end
                
                H=this.slInfo(this.AttributeRegions.urls{i});
                urls=this.reccollectUrls(H.Groups, 'AttributeValueRegion');
                
                if this.verbose>2
                    fprintf('   %d AttributeValueRegion found\n', numel(urls))
                end
                
                for j=1:numel(urls)
                    value=this.slreadatt(urls{j}, 'name');
                    if this.verbose>2
                        fprintf('     value:%s\n' , value);
                    end
                    
                    H2=this.slInfo(urls{j});
                    %seach for Regions
                    for k=1:numel(H2.Groups)
                        if regexp(H2.Groups(k).Name,'/Regions$')
                            if this.verbose>2
                                fprintf('           Region entry found\n')
                            end
                            
                            H3=this.slInfo([urls{j},'/Regions']);
                            regionidx=zeros(numel(H3.Datasets),1);
                            for l=1:numel(H3.Datasets)
                                if this.verbose>2
                                    fprintf('          -> %s\n', [H3.Name,'/', H3.Datasets(l).Name])
                                end
                                %get hexid
                                id=this.slread([H3.Name,'/', H3.Datasets(l).Name]);
                                id=id{1};
                                %seach hexid in region urls
                                ls=cellfun(@(x)~isempty(regexp(  x, [id,'$'], 'once')), this.urlRegionDataSet);
                                %be sure that result is valid 
                                idx=find(ls);
                                if isempty(idx)
                                    if this.verbose>0
                                        warning('there is no region with id <%s>\n url:%s', id, [H3.Name,'/', H3.Datasets(l).Name])
                                    end
                                    idx=0;
                                else
                                    if length(idx)>1
                                        warning('there are multiple regions with id <%s>', id);
                                    end
                                end
                                regionidx(l)=idx(1);
                                if this.verbose>2
                                    fprintf('             => region index %d\n', regionidx(l));
                                end
                            end
                                                        
                        end
                    end
                    
                    if this.verbose>2
                        fprintf('     -----------------------------------------------\n');
                        fprintf('     =>value:%s => regionindex:' , value);
                        fprintf(' %d' , regionidx);
                        fprintf('\n\n\n');
                    end
                    
                    %skip attribute values with no region set
                    regionidx=regionidx(regionidx>0);
                    
                    this.RegionProperties(regionidx, i)={value};
                    
                                       
                    
                end
                
                %type cast i-th column; empty entries are kept empty or set to nan
                switch this.AttributeRegions.ValueType(i)
                    case 1 % number
                        this.RegionProperties(:,i)=num2cell(str2double(this.RegionProperties(:,i)));
                    case 2 %binary
                        this.RegionProperties( strcmpi( this.RegionProperties(:,i), 'true') ,i)={1};
                        this.RegionProperties( strcmpi( this.RegionProperties(:,i), 'false'),i)={0};
                end
            end
        end
            
    end
    
    
     
    % low level methods
    
    methods (Access=protected)
        
              
        % recursive visit all groups with attribute 'type' == type
        % urls are returned for data acquisition with slRead
        function [urls, counter]=reccollectUrls(this, Groups, type, urls, counter)
            if nargin < 5
                counter=0;
            end
            if nargin < 4
                urls={};
            end
            
            for i=1:length(Groups)
                for j=1:length(Groups(i).Attributes)
                    key=Groups(i).Attributes(j).Name;
                    val=Groups(i).Attributes(j).Value{1};
                    %fprintf('%15s  =>  %s\n', key , val  );
                    
                    if( strcmp(key, 'type'))
                        if(strcmp(val,type))
                            counter=counter+1;
                            urls{counter}=Groups(i).Name;
                        end
                    end
                end
                
                %recursion
                if isfield(Groups(i), 'Groups')
                    [urls,counter]=this.reccollectUrls(Groups(i).Groups, type, urls, counter);
                end
                
            end
        end
        
        %return Info-struct for groupname
        %optional restriction to file scheme version
        function H=slInfo(this, groupname, default, reqver)
           if nargin < 4, reqver=-inf; end;
           if nargin < 3, default=[] ; end;
           
           if this.version >= reqver
               try
                   H=h5info(this.filename, groupname);
               catch ME
                   warning(ME.message);
                   H=default;
               end
           else
               H=default;
           end
        end
        
        
        %return dataset of given url
        %optional restriction to file scheme version
        function D=slread(this, url, default, reqver)
            if nargin < 5
                reqver=-inf;
            end
            
            if nargin < 4
                default=[];
            end
            
            if this.version >= reqver
                try
                    D=h5read(this.filename, url);
                catch ME
                    warning(ME.message)
                    D=default;
                end
            else
                D=default;
            end
        end
          
        
        
        %return data attribute of url and tag 
        %optional restriction to file scheme version
        function str=slreadatt(this, url, tag, default, reqver)
            if nargin < 5
                reqver=-inf;
            end
            
            if nargin < 4
                default='';
            end
            
            if this.version >= reqver
                try
                    str=h5readatt(this.filename, url, tag);
                    str=str{1};
                catch ME
                    warning(ME.message)
                    str=default;
                end
            else
                str=default;
            end
        end
    end
    
    
    methods(Static)
        function checkplugin
            
            if isempty(getenv('HDF5_PLUGIN_PATH')),
                
                if isunix,
                    path='../lib/Linux';
                elseif ispc,
                    path='..\lib\Windows';
                elseif ismac,
                    path='../lib/Mac';
                    warning('LZ4 plugin is not available.');
                end
                
                setenv('HDF5_PLUGIN_PATH', [fileparts(mfilename('fullpath')),'/', path]);
            end
        end
        
        function  [istart,iend] = indexset2interval(I)
            if length(I) > 1 && any(I(2:end)-I(1:end-1) ~= 1)
               error('indexset must have stride one, e.g. 1:n');
            end
           
           istart=I(1);
           iend= I(end);
           
           if istart <= 0,
               error('indices are 1-based and must be positive');
           end
           
        end
        
        function Y = indexset2intervals(I)
          % Convert index set I to list of index intervals such that
          % concatenation of all intervals is identical to I.
          % Output intervals are returned as rows of matrix Y.
          
          if isempty(I)
            Y = [];
          elseif ~(isnumeric(I) && isvector(I) && all(I>=1) && all(mod(I,1)==0))
            error('index set must be a vector of positive integers');
          else
            % Identify stretches of consecutive integers
            d=diff(I(:));
            k2=[find(d~=1); length(I)];
            k1=[1; k2(1:end-1)+1];
            Y=[I(k1) I(k2)];
          end
        end
    end
end
