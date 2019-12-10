function varargout = HoleDetector(varargin)
% HOLEDETECTOR MATLAB code for HoleDetector.fig
%      HOLEDETECTOR, by itself, creates a new HOLEDETECTOR or raises the existing
%      singleton*.
%
%      H = HOLEDETECTOR returns the handle to a new HOLEDETECTOR or the handle to
%      the existing singleton*.
%
%      HOLEDETECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HOLEDETECTOR.M with the given input arguments.
%
%      HOLEDETECTOR('Property','Value',...) creates a new HOLEDETECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HoleDetector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HoleDetector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HoleDetector

% Last Modified by GUIDE v2.5 10-Feb-2016 12:24:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HoleDetector_OpeningFcn, ...
                   'gui_OutputFcn',  @HoleDetector_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before HoleDetector is made visible.
function HoleDetector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HoleDetector (see VARARGIN)

% Choose default command line output for HoleDetector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using HoleDetector.
% if strcmp(get(hObject,'Visible'),'off')
%     plot(rand(5));
% end

% UIWAIT makes HoleDetector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HoleDetector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,filePath]=uigetfile('*.sl','Hole Detection');
ExPath = [filePath fileName];
if fileName~=0
    
    try
       handles.maldiData=MSMaldiData(ExPath);
       %handles.maldiData=handles.maldiData.reduce('single_cuts');
    catch
        set(handles.path,'String','Problems loading the .sl file');
        return
    end
    set(handles.path,'String',ExPath);
    handles.filename=ExPath;
     if ~isfield(handles,'app_started') %if it is the first time that
        handles.app_started=true(); 
        
        contents = get(handles.criterium_popupmenu,'String');
        handles.criterium=read_criterium(contents{get(handles.criterium_popupmenu, 'Value')});
        

        %making components visible
        caption = sprintf('%.2f', get(handles.radiusSlider,'Value'));
        set(handles.radiusText, 'String', caption);
        set(handles.radiusSlider,'Visible', 'on')
        set(handles.radiusText, 'Visible', 'on')

        caption = sprintf('%.2f', get(handles.thresholdSlider,'Value'));
        set(handles.thresholdText, 'String', caption);
        set(handles.thresholdSlider,'Visible', 'on')
        set(handles.thresholdText, 'Visible', 'on')

        caption = sprintf('%.2f', get(handles.radiusSliderUnder,'Value'));
        set(handles.radiusTextUnder, 'String', caption);

        caption = sprintf('%.2f', get(handles.thresholdSlider,'Value'));
        set(handles.thresholdText, 'String', caption);
        set(handles.thresholdSlider,'Visible', 'on')
        set(handles.thresholdText, 'Visible', 'on')
        
        set(handles.radiobuttonpanel, 'Visible', 'on')
        set(handles.radius_text, 'Visible', 'on')
        set(handles.threshold_text, 'Visible', 'on')
        set(handles.criterium_text, 'Visible', 'on')
        set(handles.criterium_popupmenu, 'Visible', 'on')
        set(handles.axes1, 'Visible', 'on')

        set(handles.uipanelover, 'Visible', 'on')
        set(handles.uipanelunder, 'Visible', 'on')   
        handles.holePeaks=842.5;
        handles.tissuePeaks=[];

        set(handles.uibuttongrouplogic, 'Visible', 'on')
        handles.logicFlag='&&';
        
        set(handles.checkbox_im_proc, 'Visible', 'on')
        handles.postProcFlag=get(handles.checkbox_im_proc,'Value');
     end
     guidata(hObject, handles);  
     compute_Holes(handles);
     zoom2cursor;
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)


% --- Executes on selection change in criterium_popupmenu.
function criterium_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to criterium_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns criterium_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from criterium_popupmenu
contents = get(hObject,'String');
criterium=read_criterium(contents{get(hObject, 'Value')});
if ~strcmp(criterium,handles.criterium) %if the criterium was changed recompute
    handles.criterium=criterium;
    compute_Holes(handles);
    guidata(hObject, handles);
end



% --- Executes during object creation, after setting all properties.
function criterium_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to criterium_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function radiusSlider_Callback(hObject, eventdata, handles)
% hObject    handle to radiusSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
caption = sprintf('%.2f', get(handles.radiusSlider,'Value'));
set(handles.radiusText,'String', caption);
compute_Holes(handles);


% --- Executes during object creation, after setting all properties.
function radiusSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiusSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function thresholdSlider_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
caption = sprintf('%.2f', get(handles.thresholdSlider,'Value'));
set(handles.thresholdText,'String', caption);
compute_Holes(handles);

% --- Executes during object creation, after setting all properties.
function thresholdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on radiusSlider and none of its controls.
function radiusSlider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to radiusSlider (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function compute_Holes(handles)

if ~(isempty(handles.holePeaks)&&isempty(handles.tissuePeaks))
    
    holePeakCriterium.peaks=handles.holePeaks;tissuePeakCriterium.peaks=handles.tissuePeaks;
    holePeakCriterium.radius=get(handles.radiusSlider,'Value');
    holePeakCriterium.threshold=get(handles.thresholdSlider,'Value');


    tissuePeakCriterium.radius=get(handles.radiusSliderUnder,'Value');
    tissuePeakCriterium.threshold=get(handles.thresholdSliderUnder,'Value');

    mask=MSDetectHoles(handles.maldiData, holePeakCriterium, tissuePeakCriterium, handles.logicFlag, handles.criterium, handles.postProcFlag);
    VY=0.5*(mask+1);
    I=handles.maldiData.positions.encube(VY);
    im=imagesc(I,'Parent',handles.axes1);axis off;
    set(im,'HitTest','on');
    %zoom2cursor;
    %plot(handles.axes2,handles.maldiData.mzVector,handles.maldiData.meanData), hold on
    %maldiDataHoles=handles.maldiData.reduce(mask);
    %maldiDataTissue=handles.maldiData.reduce(~mask);
    %plot(handles.axes2, maldiDataTissue.mzVector,maldiDataTissue.meanData,'r-'), hold off
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiusSlider.
function radiusSlider_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiusSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in radiobuttonpanel.
function radiobuttonpanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in radiobuttonpanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue, 'String')
    case 'raw'
        handles.maldiData.setNormalization('raw')
    case 'total ion count'
        handles.maldiData.setNormalization('tic')
    case 'median'
        handles.maldiData.setNormalization('median')
end
guidata(hObject, handles);
compute_Holes(handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
%bdfcn;
zoomparams = getappdata(gcf,'zoomparams');
% SelectionType
% normal: Click left mouse button
% extend: Shift - click left mouse button or click both left and right mouse buttons
% alt: Control - click left mouse button or click right mouse button
% open: Double click any mouse button

switch get(gcf,'selectiontype')
case 'normal'
	zoomparams.pct = max(0.01,zoomparams.pct*0.9);
case 'alt'
	zoomparams.pct = min(1,zoomparams.pct*1.1);
case 'extend'
	set(zoomparams.currax,'xlim',zoomparams.oldxlim,'ylim',zoomparams.oldylim,'zlim',zoomparams.oldzlim);
    zoomparams.pct=0.5;
case 'open'
	zoomparams.pct = max(0.01,zoomparams.pct*0.81);
end

zoomparams.xdist = zoomparams.pct*zoomparams.xrange;
zoomparams.ydist = zoomparams.pct*zoomparams.yrange;
zoomparams.zdist = zoomparams.pct*zoomparams.zrange;

setappdata(gcf,'zoomparams',zoomparams);
feval(getappdata(gcf,'zoomfcnhandle'));

return


% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)

steps=eventdata.VerticalScrollCount;
if steps
    zoomFactor=0.9^steps;
else
    zoomFactor=1.1^(-steps);
end
zoomparams = getappdata(gcf,'zoomparams');
zoomparams.pct = max(0.01,min(1,zoomparams.pct*zoomFactor));

zoomparams.xdist = zoomparams.pct*zoomparams.xrange;
zoomparams.ydist = zoomparams.pct*zoomparams.yrange;
zoomparams.zdist = zoomparams.pct*zoomparams.zrange;

setappdata(gcf,'zoomparams',zoomparams);
feval(getappdata(gcf,'zoomfcnhandle'));

return


% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
checkboxButtonDown(hObject, handles);

% --- Executes on slider movement.
function radiusSliderUnder_Callback(hObject, eventdata, handles)
% hObject    handle to radiusSliderUnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
caption = sprintf('%.2f', get(handles.radiusSliderUnder,'Value'));
set(handles.radiusTextUnder,'String', caption);
compute_Holes(handles);

% --- Executes during object creation, after setting all properties.
function radiusSliderUnder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiusSliderUnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function thresholdSliderUnder_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdSliderUnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

caption = sprintf('%.2f', get(handles.thresholdSliderUnder,'Value'));
set(handles.thresholdTextUnder,'String', caption);
compute_Holes(handles);
% --- Executes during object creation, after setting all properties.
function thresholdSliderUnder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdSliderUnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

checkboxButtonDown(hObject, handles);
% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
checkboxButtonDown(hObject, handles);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)

% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkboxButtonDown(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkboxButtonDown(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over checkbox1.
function checkbox1_ButtonDownFcn(hObject, eventdata, handles)

% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function checkboxButtonDown(hObject, handles)
chbPanelId=hObject.Parent.get('Tag');
chbMz=str2double(hObject.get('String'));
if strcmp(chbPanelId,'uipanelover')
    if hObject.get('Value')==0
        handles.holePeaks(handles.holePeaks==chbMz)=[];
    else
        handles.holePeaks=[handles.holePeaks chbMz];
    end
else
    if hObject.get('Value')==0
        handles.tissuePeaks(handles.tissuePeaks==chbMz)=[];
    else
        handles.tissuePeaks=[handles.tissuePeaks chbMz];
    end    
end
guidata(hObject, handles);
compute_Holes(handles);




% --- Executes when selected object is changed in uibuttongrouplogic.
function uibuttongrouplogic_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongrouplogic 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.logicFlag=get(eventdata.NewValue, 'String');
guidata(hObject, handles);
compute_Holes(handles);

function criterium=read_criterium(string)
switch string
    case 'norm relative to mean'
        criterium='norm_rel_mean';
    case 'sum relative to mean'
        criterium='sum_rel_mean';
    case 'norm of all features'
        criterium='norm';
    case 'sum of all features';
        criterium='sum';
end



% --- Executes on key press with focus on criterium_popupmenu and none of its controls.
function criterium_popupmenu_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to criterium_popupmenu (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_im_proc.
function checkbox_im_proc_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_im_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_im_proc
handles.postProcFlag=get(hObject,'Value');
guidata(hObject, handles);
compute_Holes(handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over checkbox_im_proc.
function checkbox_im_proc_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_im_proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15
checkboxButtonDown(hObject, handles);
