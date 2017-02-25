function varargout = visueye(varargin)
% VISUEYE MATLAB code for visueye.fig
%      VISUEYE, by itself, creates a new VISUEYE or raises the existing
%      singleton*.
%
%      H = VISUEYE returns the handle to a new VISUEYE or the handle to
%      the existing singleton*.
%
%      VISUEYE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUEYE.M with the given input arguments.
%
%      VISUEYE('Property','Value',...) creates a new VISUEYE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visueye_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visueye_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visueye

% Last Modified by GUIDE v2.5 28-Jun-2016 18:57:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visueye_OpeningFcn, ...
                   'gui_OutputFcn',  @visueye_OutputFcn, ...
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

% --- Executes just before visueye is made visible.
function visueye_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visueye (see VARARGIN)

% Choose default command line output for visueye
handles.output = hObject;
handles.good_list = [];
handles.bad_list = [];

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using visueye.
if strcmp(get(hObject,'Visible'),'off')
%     plot(rand(5));
end

% UIWAIT makes visueye wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visueye_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
% switch popup_sel_index
%     case 1
%         plot(rand(5));
%     case 2
%         plot(sin(1:0.01:25.99));
%     case 3
%         bar(1:.5:10);
%     case 4
%         plot(membrane);
%     case 5
%         surf(peaks);
% end


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

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
[handles, hObject] = cellist_LST_fun(handles, hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
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
handles = load_BTN_fun(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4


% --- Executes on button press in add_to_good.
function add_to_good_Callback(hObject, eventdata, handles)
% hObject    handle to add_to_good (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.listbox1,'String'));
selection = contents{get(handles.listbox1,'Value')};
selection = handles.df.cell_table(:,selection); % single cell
cell_struct = table2array(selection);
cell_struct.prop = handles.df.properties;
cell_struct.prop.file_name = handles.file_name;
cell_struct.cell_idx = get(handles.listbox1,'Value');
handles.good_list = [handles.good_list ;cell_struct] ;
if(any(isspace(get(handles.good_LSTBX,'String'))))
    list_contents{1} = [];
else
    list_contents{1} = cellstr(get(handles.good_LSTBX,'String'));
end
set(handles.good_LSTBX,'String',[list_contents{1};{strjoin({'Cell',num2str(cell_struct.cell_idx)},'_')}]);% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in add_to_bad.
function add_to_bad_Callback(hObject, eventdata, handles)
% hObject    handle to add_to_bad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.listbox1,'String'));
selection = contents{get(handles.listbox1,'Value')};
selection = handles.df.cell_table(:,selection); % single cell
cell_struct = table2array(selection);
cell_struct.prop = handles.df.properties;
cell_struct.prop.file_name = handles.file_name;
cell_struct.cell_idx = get(handles.listbox1,'Value');
handles.bad_list = [handles.bad_list ;cell_struct] ;
if(any(isspace(get(handles.bad_LSTBX,'String'))))
    list_contents{1} = [];
else
    list_contents{1} = cellstr(get(handles.bad_LSTBX,'String'));
end
set(handles.bad_LSTBX,'String',[list_contents{1};{strjoin({'Cell',num2str(cell_struct.cell_idx)},'_')}]);% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in good_LSTBX.
function good_LSTBX_Callback(hObject, eventdata, handles)
% hObject    handle to good_LSTBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns good_LSTBX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from good_LSTBX


% --- Executes during object creation, after setting all properties.
function good_LSTBX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to good_LSTBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bad_LSTBX.
function bad_LSTBX_Callback(hObject, eventdata, handles)
% hObject    handle to bad_LSTBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bad_LSTBX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bad_LSTBX


% --- Executes during object creation, after setting all properties.
function bad_LSTBX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bad_LSTBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stimplt_CHBX.
function stimplt_CHBX_Callback(hObject, eventdata, handles)
% hObject    handle to stimplt_CHBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stimplt_CHBX
[handles, hObject] = cellist_LST_fun(handles, handles.listbox1);
guidata(hObject, handles);

% --- Executes on button press in df_CHBX.
function df_CHBX_Callback(hObject, eventdata, handles)
% hObject    handle to df_CHBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of df_CHBX
[handles, hObject] = cellist_LST_fun(handles, handles.listbox1);
guidata(hObject, handles);

% --- Executes on button press in st_CHBX.
function st_CHBX_Callback(hObject, eventdata, handles)
% hObject    handle to st_CHBX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of st_CHBX
[handles, hObject] = cellist_LST_fun(handles, handles.listbox1);
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure();
guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles, hObject] = calc_features(handles, handles.listbox1);
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
