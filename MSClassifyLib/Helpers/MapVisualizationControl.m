function varargout = MapVisualizationControl(varargin)
% MAPVISUALIZATIONCONTROL MATLAB code for MapVisualizationControl.fig
%      MAPVISUALIZATIONCONTROL, by itself, creates a new MAPVISUALIZATIONCONTROL or raises the existing
%      singleton*.
%
%      H = MAPVISUALIZATIONCONTROL returns the handle to a new MAPVISUALIZATIONCONTROL or the handle to
%      the existing singleton*.
%
%      MAPVISUALIZATIONCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPVISUALIZATIONCONTROL.M with the given input arguments.
%
%      MAPVISUALIZATIONCONTROL('Property','Value',...) creates a new MAPVISUALIZATIONCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MapVisualizationControl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MapVisualizationControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MapVisualizationControl

% Last Modified by GUIDE v2.5 24-Nov-2016 17:19:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MapVisualizationControl_OpeningFcn, ...
                   'gui_OutputFcn',  @MapVisualizationControl_OutputFcn, ...
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

% --- Executes just before MapVisualizationControl is made visible.
function MapVisualizationControl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MapVisualizationControl (see VARARGIN)

% Choose default command line output for MapVisualizationControl
handles.output = hObject;

handles.mainAxes = axes('Parent',handles.figure1,... 
                  'Units','normalized','Position',handles.axes1.get('Position'),... 
                  'Nextplot','add','Tag','my_axes'); 
% Update handles structure
guidata(hObject, handles);

if ~isempty(varargin)&&~isa(varargin{2},'MSBasisMap')
    error('The input argument must be a valid MSBasisMap object')
else
    handles.fm = varargin{2};
end
visualString = {'heat map', 'basis vectors', 'image','basis vectors contribution (signed auc)','auc'};
if ~isfield(handles.fm.decompositionObj.addInfo,'rocMaxVals') % exclude roc options
    visualString = visualString(1:end-2);
end   
set(handles.popupmenu1, 'String', visualString);
% This sets up the initial plot - only do when we are invisible
% so window can get raised using MapVisualizationControl.
handles.axes1.set('Visible','off')
if strcmp(get(hObject,'Visible'),'off')
    axes(handles.mainAxes)
    execute(handles);
end
set(handles.nFeaturesE,'String',handles.fm.numFeatures)
sortings=get(handles.sortingPM,'String');
index=find(strcmp(handles.fm.decompositionObj.sortTypeCurr,sortings),1);
set(handles.sortingPM,'Value',index);
guidata(hObject, handles);
% UIWAIT makes MapVisualizationControl wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MapVisualizationControl_OutputFcn(hObject, eventdata, handles)
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
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

execute(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

function nFeaturesE_Callback(hObject, eventdata, handles)
% hObject    handle to nFeaturesE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nFeaturesE as text
%        str2double(get(hObject,'String')) returns contents of nFeaturesE as a double
execute(handles);

% --- Executes during object creation, after setting all properties.
function nFeaturesE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nFeaturesE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in nFeaturesPB.
function nFeaturesPB_Callback(hObject, eventdata, handles)
% hObject    handle to nFeaturesPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

execute(handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over nFeaturesPB.
function nFeaturesPB_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to nFeaturesPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
execute(handles)
function addDescription(type, handles)
string = handles.fm.visualizationDescription(type);
[outstring,~] = textwrap(handles.descriptionE,{string});
set(handles.descriptionE,'String',outstring)

function execute(handles)
axes(handles.mainAxes)
cla(handles.mainAxes)
hAxes = findobj(gcf,'type','axes');
for i=1:length(hAxes)
    cla(hAxes(i));
    axis(hAxes(i),'off');
end 
popup_sel_index = get(handles.popupmenu1,'Value');
popup_sel_string = get(handles.popupmenu1,'String');
popup_sel_string = popup_sel_string{popup_sel_index};
nFeatures = str2double(get(handles.nFeaturesE,'String'));
if isempty(nFeatures)||isnan(nFeatures)
    set(handles.nFeaturesE,'String',num2str(handles.fm.numFeatures));
    nFeatures = handles.fm.numFeatures;
end
set(handles.thresholdT,'Visible','off');
set(handles.thresholdS,'Visible','off');
modus = popup_sel_string;
switch popup_sel_string
    case 'auc'
        set(handles.thresholdT,'Visible','on');
        set(handles.thresholdS,'Visible','on');
    case 'basis vectors'
        modus='basis';
    case 'heat map'
        modus = 'heat';
    case 'basis vectors contribution (signed auc)'
        set(handles.thresholdT,'Visible','on');
        set(handles.thresholdS,'Visible','on');
        modus='signed';
end
addDescription(modus,handles);
handles.fm.show(modus,'nFeatures',nFeatures,'threshold',get(handles.thresholdS,'Value'));


% --- Executes on slider movement.
function thresholdS_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value=get(hObject,'Value');
set(handles.thresholdT,'String',['auc threshold: ' sprintf('%0.3f',value)]);
execute(handles)
    

% --- Executes during object creation, after setting all properties.
function thresholdS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in sortingPM.
function sortingPM_Callback(hObject, eventdata, handles)
% hObject    handle to sortingPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sortingPM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sortingPM
sortings=get(hObject,'String');
handles.fm.decompositionObj.sortBasis(sortings{get(hObject,'Value')});
execute(handles)

% --- Executes during object creation, after setting all properties.
function sortingPM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sortingPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function descriptionE_Callback(hObject, eventdata, handles)
% hObject    handle to descriptionE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of descriptionE as text
%        str2double(get(hObject,'String')) returns contents of descriptionE as a double


% --- Executes during object creation, after setting all properties.
function descriptionE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to descriptionE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
