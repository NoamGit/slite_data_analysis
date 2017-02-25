function varargout = visueye_2(varargin)
% VISUEYE_2 MATLAB code for visueye_2.fig
%      VISUEYE_2, by itself, creates a new VISUEYE_2 or raises the existing
%      singleton*.
%
%      H = VISUEYE_2 returns the handle to a new VISUEYE_2 or the handle to
%      the existing singleton*.
%
%      VISUEYE_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUEYE_2.M with the given input arguments.
%
%      VISUEYE_2('Property','Value',...) creates a new VISUEYE_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visueye_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visueye_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visueye_2

% Last Modified by GUIDE v2.5 25-Oct-2016 15:59:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visueye_2_OpeningFcn, ...
                   'gui_OutputFcn',  @visueye_2_OutputFcn, ...
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


% --- Executes just before visueye_2 is made visible.
function visueye_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visueye_2 (see VARARGIN)

% Choose default command line output for visueye_2
set(handles.param_SLDR, 'min', 0);
set(handles.param_SLDR, 'max', 10);
set(handles.param_SLDR, 'SliderStep', [0.001,0.01]);
set(handles.param_SLDR, 'Value', 0); % Somewhere between max and min.
set(handles.itr_BTN,'enable','off');
set(handles.thresvalue_TXT,'String',num2str(get(handles.param_SLDR,'Value')));
handles.output = hObject;
handles.nographics = 0;
handles.dest_path = [];
set(handles.thresvalue_TXT,'String',num2str(get(handles.param_SLDR,'Value')));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visueye_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visueye_2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in cell_location_LST.
function cell_location_LST_Callback(hObject, eventdata, handles)
% hObject    handle to cell_location_LST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell_location_LST contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell_location_LST
[handles, hObject] = list_short_fun(handles, handles.cell_location_LST);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cell_location_LST_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_location_LST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cell_LST.
function cell_LST_Callback(hObject, eventdata, handles)
% hObject    handle to cell_LST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell_LST contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell_LST
mode = get(handles.auto_CBX,'Value');
[handles, hObject] = list_long_fun(handles, handles.cell_LST, mode);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cell_LST_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_LST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_BTN.
function load_BTN_Callback(hObject, eventdata, handles)
% hObject    handle to load_BTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = load_fn(handles,hObject);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in F_CBX.
function F_CBX_Callback(hObject, eventdata, handles)
% hObject    handle to F_CBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F_CBX
mode = 0;
[handles, hObject] = list_long_fun(handles, handles.cell_LST, mode);
guidata(hObject, handles);

% --- Executes on button press in path_BTN.
function path_BTN_Callback(hObject, eventdata, handles)
% hObject    handle to path_BTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dest_path = uigetdir('C:\Users\noambox\Documents\Sync\Neural data\stim_features','Save to..');
guidata(hObject, handles);

% --- Executes on button press in C_CBX.
function C_CBX_Callback(hObject, eventdata, handles)
% hObject    handle to C_CBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of C_CBX
mode = 0;
[handles, hObject] = list_long_fun(handles, handles.cell_LST, mode);
guidata(hObject, handles);

% --- Executes on button press in S_CBX.
function S_CBX_Callback(hObject, eventdata, handles)
% hObject    handle to S_CBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of S_CBX
mode = 0;
[handles, hObject] = list_long_fun(handles, handles.cell_LST, mode);
guidata(hObject, handles);

% --- Executes on button press in save_BTN.
function save_BTN_Callback(hObject, eventdata, handles)
% hObject    handle to save_BTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hold all;
axes(handles.axes2);
plot(randi(10,100,1));
guidata(hObject, handles);


function fName_ETXT_Callback(hObject, eventdata, handles)
% hObject    handle to fName_ETXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fName_ETXT as text
%        str2double(get(hObject,'String')) returns contents of fName_ETXT as a double


% --- Executes during object creation, after setting all properties.
function fName_ETXT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fName_ETXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=figure;ax=axes;clf;
new_handle = copyobj(handles.axes2,fig); % Copy axes object h into figure f1
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=figure;ax=axes;clf;
new_handle = copyobj(handles.axes1,fig); % Copy axes object h into figure f1
% set(gca,'ActivePositionProperty','outerposition')
% set(gca,'Units','normalized')
% set(gca,'OuterPosition',[0 0 1 1])
set(gcf,'position',[2100 0 871 250])
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=figure;ax=axes;clf;
new_handle = copyobj(handles.axes4,fig); % Copy axes object h into figure f1
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5


% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=figure;ax=axes;clf;
new_handle = copyobj(handles.axes5,fig); % Copy axes object h into figure f1
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])
guidata(hObject, handles);


% --- Executes on slider movement.
function param_SLDR_Callback(hObject, eventdata, handles)
% hObject    handle to param_SLDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.thresh_TXT,'String',['Spike thresh value: ',num2str(get(handles.param_SLDR,'Value'))]);
mode = handles.auto_CBX.Value;
[handles, hObject] = list_long_fun(handles, handles.cell_LST, mode);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function param_SLDR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param_SLDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in eval_BTN.
function eval_BTN_Callback(hObject, eventdata, handles)
% hObject    handle to eval_BTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mode = 1;
[handles, hObject] = list_long_fun(handles, handles.cell_LST, mode);
guidata(hObject, handles);


function param_ETXT_Callback(hObject, eventdata, handles)
% hObject    handle to param_ETXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of param_ETXT as text
%        str2double(get(hObject,'String')) returns contents of param_ETXT as a double


% --- Executes during object creation, after setting all properties.
function param_ETXT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param_ETXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_CBX.
function auto_CBX_Callback(hObject, eventdata, handles)
% hObject    handle to auto_CBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_CBX


% --------------------------------------------------------------------
function uipanel8axes2_parent_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel8axes2_parent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in athresh_CHBX.
function athresh_CHBX_Callback(hObject, eventdata, handles)
% hObject    handle to athresh_CHBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of athresh_CHBX


% --- Executes on button press in plotmore_CHBX.
function plotmore_CHBX_Callback(hObject, eventdata, handles)
% hObject    handle to plotmore_CHBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotmore_CHBX
% --- Executes on button press in athresh_CHBX.
function scale_CHBX_Callback(hObject, eventdata, handles)
% hObject    handle to athresh_CHBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of athresh_CHBX


% --- Executes on button press in loadF_BTN.
function loadF_BTN_Callback(hObject, eventdata, handles)
% hObject    handle to loadF_BTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
funcName = uigetfile(['C:\Users\noambox\Documents\Sync\code\Data Analysis\functions\*.m'],'Select function');
handles.funcName = funcName(1:end-2); 
set(handles.itr_BTN,'enable','on');
guidata(hObject, handles);

% --- Executes on button press in itr_BTN.
function itr_BTN_Callback(hObject, eventdata, handles)
% hObject    handle to itr_BTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = callIterator_fun(handles);
if ~isempty(handles.dest_path)   
    table2save = handles.feature_table;
    save([handles.dest_path,'\crp_last_analysis_clust.mat'],'table2save');
end
guidata(hObject, handles);
