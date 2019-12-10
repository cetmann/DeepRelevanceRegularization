function MSGenerateReport(title,maldiData,varargin)

% Generate a report for an MSMaldiData object.
%
% MSGenerateReport(title, maldiData, varargin):
% 
% Mandatory input arguments:
%
% title - A string containing the title of the document/data set.
%
% maldiData - The MSMaldiData object for which the report will be created.
%
% Optional input arguments:
%
% 'draft' - Create a rough draft of your report. This option makes report
% generation faster.
%
% 'skipQuality' - skip calculation of the quality measure and pictures
%
% 'FontName', string - A string containing the desired fonts name. If the input does
% not match any of the system fonts, the function tries to find the right name.
% If this fails, Arial is used instead. Best results are archieved by supplying
% the correct system font name right away.
%
% 'FontSize', int/string - Specify the font size for the report. To specify a
% font size of 16pt, input could be: 16, '16' or '16pt'
%
% 'savePath', string - Specify the save path for your file 


%% DEFAULT VALUES
fontName = 'Myriad Pro Light';
fontSize = '12pt';
draft = false;
skipQuality = false;
path = '';

qualityMeasure = 1.0123456789;
sNR = pi;


%% PROCESS AND VALIDATE INPUT
if ~ischar(title)
    throw(MException('MSClassifyLib:MSGenerateReport:InvalidTitleInput',...
        'The title string input for MSGenerateReport() must be of type char.'));
end

if ~isa(maldiData, 'MSData')
    throw(MException('MSClassifyLib:MSGenerateReport:InvalidMaldiDataInput',...
        'The data input for MSGenerateReport() must be an MSData object.'));
end

mandatoryArguments = -nargin('MSGenerateReport')-1;
if nargin > mandatoryArguments
    iArgument = 1;
    while iArgument <= nargin - mandatoryArguments
        currArg = varargin{iArgument};
        if isempty(currArg)
            currArg = '';
        end
        
        switch lower(currArg)
            case 'fontname'
                fontName = varargin{iArgument+1};
                iArgument = iArgument + 2;
            case 'fontsize'
                fS = varargin{iArgument+1};
                if isnumeric(fS)
                    fontSize = [num2str(fS) 'pt'];
                else
                    if ischar(fS)
                        if strcmpi('pt',fS(end-1:end))
                            fS = fS(1:end-2);
                        end
                        if ~isempty(str2num(fS))
                            fontSize = [fS 'pt'];
                        else
                            warning('Your input for the variable fontSize seems to be invalid. A default value is used instead.');
                        end
                    end
                end
                fS = str2num(fontSize(1:end-2));
                if fS < 8 || fS > 25
                    warning('Your input for the variable fontSize seems unreasonable. Acceptable values lie in the interval 8pt - 25pt. The default font size is used instead.');
                end
                
                iArgument = iArgument + 2;
            case 'draft'
                draft = true;
                iArgument = iArgument + 1;
            case 'skipquality'
                skipQuality = true;
                iArgument = iArgument + 1;
            case 'savepath'
                input = varargin{iArgument + 1};
                if ischar(input) && exist(input,'dir') == 7
                    %Append missing '/'
                    if ~strcmp(input(end),'/')
                       input(end + 1)  = '/';
                    end
                    path = input;
                    iArgument = iArgument + 2;
                else
                    throw(MException('MSClassifyLib:MSGenerateReport:InvalidSavePath',...
                    'The input save path must be of type char and a valid directory on your machine.'));
                end
            otherwise
                warning(['Unknown argument: ' currArg]);
                iArgument = iArgument + 1;
        end
    end
end

%% SET-UP
%Initialize number to keep track of plots
plotNumber = 0;

% Import report generator and create document
import mlreportgen.dom.*;

if ~isempty(path)
    d = Document([path title 'Report'],'docx');
else
    d = Document([title 'Report'],'docx');
end


%Suppress plot output from now on
fig = figure('Visible','off');
set(gca,'position',[.05 .05 .9 .9],'units','normalized');

%Calculate variables
if ~draft
    %Calculate Signal-To-Noise-Ratio of mean spectrum
    fprintf(1,'Calculating SNR \t\t…\n');
    sNR = max(maldiData.meanData)/median(maldiData.meanData);
    fprintf(1,'\b\bDone\n');
    
    %Generate max picture
    fprintf(1,'Calculating maximum picture \t…\n');
    maxPic = zeros(numel(maldiData.positions.indexGrid),1);
    maxPic(logical(maldiData.positions.indexGrid),1) = max(maldiData.data,[],2);
    maxPic = reshape(maxPic,size(maldiData.positions.indexGrid,1),...
        size(maldiData.positions.indexGrid,2),[]);
    fprintf(1,'\b\bDone\n');
    
    %Calculating quality measure
    if ~skipQuality
        fprintf(1,'Calculating quality measure \t…\n');
        %imStats = getImStats(maldiData.data,maldiData.positions.indexGrid);
        imStats = evalin('base','imStats');
        
        %Calculate the quality measure and stats. Suppress console output.
        [~,qualityMeasure, qStats] = evalc('getQuality(maldiData.data,imStats);');
        
        fprintf(1,'Calculating quality measure \tDone\n');
    else
        fprintf(1,'Skipping quality measure calculation\n');
    end
end

fprintf(1,'Gathering report resources \t…\n');
%Get current MATLAB version. This is used later, because some features only work
%with R2015b and newer versions.
charCode = @(c) find(char(0:255) == c) - 1;
matlabVersion = ver;
release = matlabVersion.Release;
releaseYear = str2double(release(3:end-2));
releaseLetter = charCode(release(end-1));
newerThanR2015b = (releaseYear == 2015 && releaseLetter >= 98) || releaseYear > 2015;


% Find out, if a font with the desired name is installed. If so, use this font
% for the Report. Else, use the standard font 'Arial'.
fontsList = listfonts;
if any(cell2mat(strfind(fontsList,fontName)));
    fontIndices = strfind(fontsList,fontName);
    for k = 1:size(fontIndices,1)
        if isempty(fontIndices{k,1})
            fontIndices{k,1} = 0;
        end
    end
    fontIndices = find([fontIndices{:}]==1);
    names = fontsList(fontIndices);
    lengths=cellfun(@(x) numel(x),names);
    systemFontName=names(lengths==min(lengths));
    systemFontName = systemFontName{:};
else
    systemFontName = 'Arial';
end

%Define Styles
global titleStyle dateAndVersionStyle subtitleStyle textStyle numberStyle
titleStyle = {Color('Gray'), FontFamily(systemFontName),...
    FontSize(changedFontSize(13)), HAlign('center')};
dateAndVersionStyle = {Color('LightGray'), FontFamily(systemFontName),...
    FontSize(changedFontSize(-4)), HAlign('right')};
subtitleStyle = {Color('Gray'),FontFamily(systemFontName),...
    FontSize(changedFontSize(4))};
textStyle = {Color('Gray'), FontFamily(systemFontName), FontSize(fontSize)};
numberStyle = {FontFamily(systemFontName), FontSize(fontSize), Italic(true)};
fprintf(1,'\b\bDone\n');
%% GENERATE REPORT
fprintf(1,'Generating report \t\t…\n');

% Title
appendPar(title,titleStyle);

% Creation date
creationDate = datetime('today');
creationDate = datestr(creationDate,29);
creationDatePar = Paragraph(['created ', creationDate]);
creationDatePar.HAlign = 'right';
creationDatePar.Color = 'lightGray';
append(creationDatePar, [' using MATLAB version ' release(2:end-1)]);
appendPar(creationDatePar,dateAndVersionStyle);

%General information
appendPar('General information',subtitleStyle);

appendNumberPar('\tNumber of spectra: \t', num2str(maldiData.numItems),...
    '\t\tNumber of m/z bins: \t', num2str(maldiData.dataLength));

appendNumberPar('\tM/z resolution: \t', num2str(maldiData.mzResolution),...
    '\t\tNormalization: \t', num2str(maldiData.normalization));
newLine();

% Mean spectrum
appendPar('Mean spectrum',subtitleStyle);
meanSpecPlot = plot(maldiData.mzVector,maldiData.meanData);
axis tight;
appendPlot(meanSpecPlot,6,2,'center');
newLine();

% Maximum picture
appendPar('Maximum picture',subtitleStyle);
if draft
    maxPicPlot = imshow(imread('tissue.png'));
else
    maxPicPlot = imagesc(maxPic,'AlphaData',logical(maldiData.positions.indexGrid),...
        prctile(maxPic(:), [0 99.5]));
    axis off;
end
appendPlot(maxPicPlot,4,3,'center');
newLine();

%Quality measures
appendPar('Quality',subtitleStyle);
appendNumberPar('\tSignal-To-Noise-Ratio: \t\t', num2str(sNR), '\t\tQuality measure: \t', num2str(qualityMeasure,10));
newLine();

% Quality pictures
appendPar('MOC-histogram and sigmoidal',subtitleStyle);
if ~draft && ~skipQuality
    MOCPic = plotMOCHistAndSigmoidal(imStats,qStats);
else
    MOCPic = imshow(imread('tissue.png'));
end
appendPlot(MOCPic,4,3,'center');
newLine();

appendPar('Sample maximum correlation, histogram and w-function',subtitleStyle);
if ~draft && ~skipQuality
    wFunPic = plotWFunAndMaxCorr(qStats);
else
    wFunPic = imshow(imread('tissue.png'));
end
appendPlot(wFunPic,4,3,'center');

if close(d);
    fprintf(1,'\b\bDone\n');
end
%% EXPORT REPORT
% Export as a DOCX-file on Unix and as a PDF on Windows
if ispc
    rptview(d.OutputPath,'pdf');
    fprintf(['\nYour report has been saved to the location:\n', d.OutputPath, '\n']);
else
    fprintf('\nYour platform does not support automatic PDF generation. Your report is generated as a DOCX-file instead. If you wish output as a PDF-file, please change the format of the report manually.\n');
    fprintf(['\nYour report has been saved to the location:\n', d.OutputPath, '\n\n']);
end

%% REMOVE UNNECESSARY FILES AND CLOSE FIGURE
for k = 1 : plotNumber
    fileName = ['reportPlot' num2str(k) '.png'];
    if exist(fileName,'file') == 2
        delete(fileName);
    end
end

close(fig);
%% AUXILIARY FUNCTIONS
    function appendPar(par, style)
        import mlreportgen.dom.*;
        
        if ~isa(par,'mlreportgen.dom.Paragraph') && ischar(par)
            par = Paragraph(par);
        end
        if exist('style','var')
            par.Style = style;
        end
        append(d,par);
    end

    function appendNumberPar(text1, italic1, text2, italic2)
        import mlreportgen.dom.*;
        par = Paragraph();
        
        text = Text(sprintf(text1));
        text.Style = textStyle;
        append(par,text);
        
        number = Text(italic1);
        number.Style = numberStyle;
        append(par,number);
        
        text = Text(sprintf(text2));
        text.Style = textStyle;
        append(par,text);
        
        number = Text(italic2);
        number.Style = numberStyle;
        append(par,number);
        
        append(d,par);
    end

    function plotHorizontalLine()
        
        %Plot a horizontal line, if the current matlab version is R2015b or newer
        if newerThanR2015b
            hr = HorizontalRule();
            append(d,hr);
        end
    end

    function newLine()
        import mlreportgen.dom.*;
        append(d,Paragraph(''));
    end

    function appendPlot(plotName, width, height, hAlign)
        import mlreportgen.dom.*;
        plotNumber = plotNumber + 1;
        
        name = ['reportPlot' num2str(plotNumber)];
        imgName= [name '.png'];
        
        saveas(plotName,imgName);
        img = Image(imgName);
        if ~exist('width','var')
            img.Width = '5in';
        else
            img.Width = [num2str(width) 'in'];
        end
        if ~exist('height','var')
            img.Height = '4in';
        else
            img.Height = [num2str(height) 'in'];
        end
        
        par = Paragraph(img);
        if exist('hAlign','var')
            par.HAlign = hAlign;
        end
        append(d,par);
    end

    function newFontSize = changedFontSize(change)
        newFontSize = [num2str(str2num(fontSize(1:end-2)) + change) 'pt'];
    end

    function p = plotMOCHistAndSigmoidal(imStats,qStats)
        [distrMoc.bars,distrMoc.centers] = hist(imStats.mocVals,30);
        distrMoc.rescaledBars = distrMoc.bars/max(distrMoc.bars);
        bar(distrMoc.centers,distrMoc.rescaledBars,0.5);
        hold on;
        xx1 = linspace(0,qStats.mocThr,300);
        xx2 = linspace(qStats.mocThr,1.05*max(distrMoc.centers),300);
        yy1 = qStats.sigmoidal(xx1);
        yy2 = qStats.sigmoidal(xx2);
        p = plot(xx1,yy1,'r','LineWidth',2);
        plot(xx2,yy2,'r:','LineWidth',2);
        axis tight;
        line([qStats.mocThr qStats.mocThr],[0 1],...
            'LineWidth',3,'Color',[1 0 0]);
        %title('Moc-histogram and sigmoidal');
        axP = gca;
        axP.FontSize = 10;
        axP.FontWeight = 'bold';
    end

    function p = plotWFunAndMaxCorr(qStats)
        if ~isempty(qStats.corrExBarC) && ~isempty(qStats.corrExBarH)
            bar(qStats.corrExBarC,qStats.corrExBarH,'FaceColor',[0.9 0.9 0.9],...
                'BarWidth',0.8);
        end
        hold on;
        
        xx = -0.2:.01:1;
        yy = qStats.wFun(xx);
        
        p = plot(xx,yy,'b','LineWidth',2);
        axP = gca;
        axP.FontSize = 10;
        axP.FontWeight = 'bold';
        axis tight;


        hxL = xlabel('Correlation');
        hxL.FontSize = 12;
        hxL.FontWeight = 'bold';
        hyL = ylabel('weight w(Corr.)');
        hyL.FontSize = 12;
        hyL.FontWeight = 'bold';
                hold off
    end

end