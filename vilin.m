function varargout = vilin(varargin)
% VILIN MATLAB code for vilin.fig
%      VILIN, by itself, creates a new VILIN or raises the existing
%      singleton*.
%
%      H = VILIN returns the handle to a new VILIN or the handle to
%      the existing singleton*.
%
%      VILIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VILIN.M with the given input arguments.
%
%      VILIN('Property','Value',...) creates a new VILIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the VILIN before vilin_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vilin_OpeningFcn via varargin.
%
%      *See VILIN Options on GUIDE's Tools menu.  Choose "VILIN allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vilin

% Last Modified by GUIDE v2.5 25-Jun-2017 14:41:41

% Begin initialization code - DO NOT EDIT

if ispc
  addpath(strcat('GUI',filesep,'Windows'));
else
  addpath(strcat('GUI',filesep,'Linux'));
end

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vilin_OpeningFcn, ...
                   'gui_OutputFcn',  @vilin_OutputFcn, ...
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


% --- Executes just before vilin is made visible.
function vilin_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vilin (see VARARGIN)

% Choose default command line output for vilin

% Custom classes added to path

 set(hObject,'Resize','off'); % disable resizing

addpath Util
% 'Functions' folder added to path
addpath(genpath('Functions'));
% Methods added to path
addpath(genpath('Methods'));
handles.presets = presets();
handles.GuiHelpers = GuiHelpers();
guidata(hObject,handles)
%adding title and labels for graphics
handles.GuiHelpers.setGraphicTitles(handles);
handles.enabledLineSearch = 'All';

INITIAL_DIMENSION = 100;

%-------adding default starting points-----------
handles.startingPoints = StartingPointGenerator(INITIAL_DIMENSION);
guidata(hObject, handles);

handles.output = hObject;

%Reading onedimensional functions
functions=dir('Functions/OneDimensional');
set(handles.oneDimFunctionPopUp, 'String', handles.GuiHelpers.fNames(functions));

%Reading multidimensional functions
functions=dir('Functions/MultiDimensional');
set(handles.multFunctionPopUp, 'String',handles.GuiHelpers.fNames(functions));

%Reading onedimensional methods
oneMethods=dir('Methods/OneDimensional');
set(handles.oneDimMethodPopUp, 'String', handles.GuiHelpers.fNames(oneMethods));

%Reading multidimensional methods
multMethods=dir('Methods/MultiDimensional');
set(handles.multiDimMethodPopUp, 'String', handles.GuiHelpers.fNames(multMethods));

%Reading method groups
methodGroups=dir('Methods/MultiDimensional');
methodGroups = methodGroups(strcmp({methodGroups.name},'LineSearch') == 0);
set(handles.methodGroupPopUp, 'String', handles.GuiHelpers.subDirs(methodGroups));
%set(handles.methodNamePopUp, 'String', subDirs(methodGroups));

%Reading line search methods
lineSearchMethods=dir('Methods/MultiDimensional/LineSearch');%Line Search methods are not standalone
set(handles.lineSearchPopUp, 'String', handles.GuiHelpers.fNames(lineSearchMethods));

%Hide calculating panel
set(handles.calculatingPanel, 'Visible', 'Off');

% Initialy set line search method
handles.lineSearchMethod = handles.GuiHelpers.getCurrentPopupString(handles.lineSearchPopUp);

% Update handles structure

guidata(hObject, handles);

%initially set starting point
multFunctionPopUp_Callback(hObject, eventdata, handles);

%initially set method group
methodGroupPopUp_Callback(hObject, eventdata, handles);

%initially activate onedimensional part of app
multiDim_Callback(hObject, eventdata, handles);

% UIWAIT makes vilin wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Executes when user attempts to close vilin.
function vilin_CloseRequestFcn(hObject, eventData, handles)
% hObject    handle to vilin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
if ispc
  rmpath(strcat('GUI',filesep,'Windows'));
else
  rmpath(strcat('GUI',filesep,'Linux'));
end

% --- Outputs from this function are returned to the command line.
function varargout = vilin_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in multiDimMethodPopUp.
function multiDimMethodPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to multiDimMethodPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns multiDimMethodPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from multiDimMethodPopUp
handles.GuiHelpers.crtaj_panel(handles)
method = handles.GuiHelpers.getCurrentPopupString(handles.multiDimMethodPopUp);
methodGroup = handles.GuiHelpers.getCurrentPopupString(handles.methodGroupPopUp);
defaultLineSearchPos = handles.GuiHelpers.getDefaultLineSearchPos(handles, method, methodGroup);
set(handles.lineSearchPopUp, 'Value', defaultLineSearchPos);
handles.lineSearchMethod = handles.GuiHelpers.getCurrentPopupString(handles.lineSearchPopUp);
handles.GuiHelpers.setLineSearchParams(handles, handles.lineSearchMethod);
guidata(hObject, handles);
handles.enabledLineSearch = handles.GuiHelpers.enableLineSearch(handles, method);
guidata(hObject, handles);

% if strcmp(handles.GuiHelpers.enableLineSearch(handles, method), 'None') == 1 || ...
%    strcmp(handles.GuiHelpers.enableLineSearch(handles, methodGroup), 'None') == 1
%     set(handles.lineSearchPopUp, 'Visible', 'Off');
% else
%     set(handles.lineSearchPopUp, 'Visible', 'On');
% end

if ~handles.GuiHelpers.enableAdvancedPanel(handles, methodGroup) || ...
   ~handles.GuiHelpers.enableAdvancedPanel(handles, method)
   set(handles.cetiriPromenljive_panel, 'Visible', 'Off');
else
   set(handles.cetiriPromenljive_panel, 'Visible', 'On');
end

set(handles.stop_cond_panel, 'Visible', 'On');


% --- Executes during object creation, after setting all properties.
function multiDimMethodPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to multiDimMethodPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function startingPointEdit_Callback(hObject, eventdata, handles)
% hObject    handle to startingPointEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startingPointEdit as text
%        str2double(get(hObject,'String')) returns contents of startingPointEdit as a double
newSP = str2num(eventdata.Source.String);
set(handles.dimensionEdit, 'String', length(newSP));


% --- Executes during object creation, after setting all properties.
function startingPointEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startingPointEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stepSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to stepSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of stepSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function stepSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dimensionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dimensionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dimensionEdit as text
%        str2double(get(hObject,'String')) returns contents of dimensionEdit as a double
newDim = str2num(eventdata.Source.String);
handles.GuiHelpers.updateStartingPoint(handles, handles.GuiHelpers.selectedFunction(handles.multFunctionPopUp), newDim);

% --- Executes during object creation, after setting all properties.
function dimensionEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dimensionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stepNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to stepNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepNumEdit as text
%        str2double(get(hObject,'String')) returns contents of stepNumEdit as a double


% --- Executes during object creation, after setting all properties.
function stepNumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to epsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsEdit as text
%        str2double(get(hObject,'String')) returns contents of epsEdit as a double


% --- Executes during object creation, after setting all properties.
function epsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in multFunctionPopUp.
function multFunctionPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to multFunctionPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns multFunctionPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from multFunctionPopUp
handles.GuiHelpers.prepParams(handles);

% --- Executes during object creation, after setting all properties.
function multFunctionPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to multFunctionPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fminResultEdit_Callback(hObject, eventdata, handles)
% hObject    handle to fminResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fminResultEdit as text
%        str2double(get(hObject,'String')) returns contents of fminResultEdit as a double


% --- Executes during object creation, after setting all properties.
function fminResultEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fminResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function iterNumResultEdit_Callback(hObject, eventdata, handles)
% hObject    handle to iterNumResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterNumResultEdit as text
%        str2double(get(hObject,'String')) returns contents of iterNumResultEdit as a double


% --- Executes during object creation, after setting all properties.
function iterNumResultEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterNumResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in multDimMinimize.
function multDimMinimize_Callback(hObject, eventdata, handles)
% hObject    handle to multDimMinimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear previous results
handles.GuiHelpers.clearMultDimResults(handles);
handles.GuiHelpers.clearGraphics(handles);
%load function and method
SELECTED_FUNCTION = handles.GuiHelpers.selectedFunction(handles.multFunctionPopUp);
SELECTED_METHOD = handles.GuiHelpers.selectedMethod(handles.multiDimMethodPopUp);
%load method params
[methodParams, success, msg] = handles.GuiHelpers.loadMultDimMethodParams(handles, SELECTED_FUNCTION);
if ~success
    msgbox(msg);
    return;
end
%show notification
set(handles.calculatingPanel, 'Visible', 'On', 'Position', [895, 460, 380, 30]);
pause(0.0001);

%calculating
try
[ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = SELECTED_METHOD(SELECTED_FUNCTION, methodParams);
catch ex
    set(handles.calculatingPanel, 'Visible', 'Off');
    
    if(strcmp(ex.identifier, 'NumOpt:implementationError'))
        msgbox(ex.message, 'Error');
        return;
    else
        msgbox('An error occured, check console output for more information.', 'Error');
        rethrow(ex);
    end
end

results = Results(fmin, xmin, valuesPerIter.gradientPerIteration(end), iterNum, cpuTime, evalNumbers, valuesPerIter);
%display and plot results
handles.GuiHelpers.displayMultDimResults(results, handles);
handles.GuiHelpers.plotMultDimResults(results, handles);
%=======================
%needed for ploting on button click
handles.iterations = valuesPerIter.iterations;
handles.gradPerIter = valuesPerIter.gradientPerIteration;
handles.valuesPerIter = valuesPerIter.functionPerIteration;
guidata(hObject, handles);
%=======================
%set values for later plotting gradient
handles.GuiHelpers.initSliders(handles, length(valuesPerIter.iterations));
%remove notification
set(handles.calculatingPanel, 'Visible', 'Off');


% --- Executes on button press in gradPlot.
function gradPlot_Callback(hObject, eventdata, handles)
% hObject    handle to gradPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.GuiHelpers.plotGrad(handles);


function cpuTimeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to cpuTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpuTimeEdit as text
%        str2double(get(hObject,'String')) returns contents of cpuTimeEdit as a double


% --- Executes during object creation, after setting all properties.
function cpuTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpuTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function funEvalNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to funEvalNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of funEvalNumEdit as text
%        str2double(get(hObject,'String')) returns contents of funEvalNumEdit as a double


% --- Executes during object creation, after setting all properties.
function funEvalNumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to funEvalNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gradEvalNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gradEvalNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gradEvalNumEdit as text
%        str2double(get(hObject,'String')) returns contents of gradEvalNumEdit as a double


% --- Executes during object creation, after setting all properties.
function gradEvalNumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gradEvalNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gradPlotStart_Callback(hObject, eventdata, handles)
% hObject    handle to gradPlotStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gradPlotStart as text
%        str2double(get(hObject,'String')) returns contents of gradPlotStart as a double
v = handles.GuiHelpers.getVector(handles.gradPlotStart);
minVal = get(handles.gradStartSlider, 'Min');
maxVal = get(handles.gradStartSlider, 'Max');
if(v >= minVal && v <= maxVal)
    set(handles.gradStartSlider, 'Value', v);
    gradStartSlider_Callback(hObject, eventdata, handles);
end


% --- Executes during object creation, after setting all properties.
function gradPlotStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gradPlotStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gradPlotEnd_Callback(hObject, eventdata, handles)
% hObject    handle to gradPlotEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gradPlotEnd as text
%        str2double(get(hObject,'String')) returns contents of gradPlotEnd as a double
v = handles.GuiHelpers.getVector(handles.gradPlotEnd);
minVal = get(handles.gradEndSlider, 'Min');
maxVal = get(handles.gradEndSlider, 'Max');
if(v >= minVal && v <= maxVal)
    set(handles.gradEndSlider, 'Value', v);
    gradEndSlider_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function gradPlotEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gradPlotEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function gradStartSlider_Callback(hObject, eventdata, handles)
% hObject    handle to gradStartSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = round(get(handles.gradStartSlider, 'value'));
set(handles.gradPlotStart, 'Value', v);
set(handles.gradPlotStart, 'String', v);
endSliderMax = get(handles.gradEndSlider, 'Max');
set(handles.gradEndSlider, 'Min', min(v+1, endSliderMax-1));
if round(get(handles.gradEndSlider, 'value')) < v + 1
    set(handles.gradEndSlider, 'Value', min(v + 1, endSliderMax));
end

gradPlot_Callback(hObject, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function gradStartSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gradStartSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function gradEndSlider_Callback(hObject, eventdata, handles)
% hObject    handle to gradEndSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = round(get(handles.gradEndSlider, 'value'));
set(handles.gradPlotEnd, 'Value', v);
set(handles.gradPlotEnd, 'String', v);
set(handles.gradStartSlider, 'Max', v);
set(handles.gradStartSlider, 'Max', max(2, v-1));
if round(get(handles.gradStartSlider, 'value')) > v - 1
    set(handles.gradStartSlider, 'Value', max(1, v - 1));
end

gradPlot_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function gradEndSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gradEndSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function deltaMin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to deltaMin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltaMin_edit as text
%        str2double(get(hObject,'String')) returns contents of deltaMin_edit as a double


% --- Executes during object creation, after setting all properties.
function deltaMin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaMin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function od_edit_Callback(hObject, eventdata, handles)
% hObject    handle to od_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of od_edit as text
%        str2double(get(hObject,'String')) returns contents of od_edit as a double


% --- Executes during object creation, after setting all properties.
function od_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to od_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function do_edit_Callback(hObject, eventdata, handles)
% hObject    handle to do_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of do_edit as text
%        str2double(get(hObject,'String')) returns contents of do_edit as a double


% --- Executes during object creation, after setting all properties.
function do_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to do_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in advanceParameter_checkbox.
function advanceParameter_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to advanceParameter_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of advanceParameter_checkbox
if get(hObject,'Value')
        set(handles.lineSearch_label, 'Visible', 'on');
        set(handles.lineSearchPopUp, 'Visible', 'on');
        set(handles.beta_edit, 'Visible','on');
        set(handles.beta_text, 'Visible','on');
        set(handles.sigma_edit, 'Visible','on');
        set(handles.sigma_text, 'Visible','on');
        set(handles.rho_edit, 'Visible','on');
        set(handles.rho_text, 'Visible','on');
        set(handles.startingPoint_edit, 'Visible','on');
        set(handles.startingPoint_text, 'Visible','on');
        %set(handles.ksi_edit, 'Visible','on');
        %set(handles.ksi_text, 'Visible','on');
        set(handles.M_edit, 'Visible','on');
        set(handles.M_text, 'Visible','on');
else
        set(handles.lineSearch_label, 'Visible', 'off');
        set(handles.lineSearchPopUp, 'Visible', 'off');
        set(handles.beta_edit, 'Visible','off');
        set(handles.beta_text, 'Visible','off');
        set(handles.sigma_edit, 'Visible','off');
        set(handles.sigma_text, 'Visible','off');
        set(handles.rho_edit, 'Visible','off');
        set(handles.rho_text, 'Visible','off');
        set(handles.startingPoint_edit, 'Visible','off');
        set(handles.startingPoint_text, 'Visible','off');    
        set(handles.ksi_edit, 'Visible','off');
        set(handles.ksi_text, 'Visible','off');
        set(handles.M_edit, 'Visible','off');
        set(handles.M_text, 'Visible','off');
end


function beta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to beta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_edit as text
%        str2double(get(hObject,'String')) returns contents of beta_edit as a double


% --- Executes during object creation, after setting all properties.
function beta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_edit as text
%        str2double(get(hObject,'String')) returns contents of sigma_edit as a double


% --- Executes during object creation, after setting all properties.
function sigma_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rho_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rho_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rho_edit as text
%        str2double(get(hObject,'String')) returns contents of rho_edit as a double


% --- Executes during object creation, after setting all properties.
function rho_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rho_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function startingPoint_edit_Callback(hObject, eventdata, handles)
% hObject    handle to startingPoint_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startingPoint_edit as text
%        str2double(get(hObject,'String')) returns contents of startingPoint_edit as a double


% --- Executes during object creation, after setting all properties.
function startingPoint_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startingPoint_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in stop_cond_checkbox.
function stop_cond_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to stop_cond_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop_cond_checkbox
if get(hObject,'Value')
        set(handles.stepNum_edit, 'Visible','on');
        set(handles.stepNum_label, 'Visible','on');
        set(handles.eps_edit, 'Visible','on');
        set(handles.eps_label, 'Visible','on');
        set(handles.workPrec_edit, 'Visible','on');
        set(handles.workPrec_label, 'Visible','on');
else
        set(handles.stepNum_edit, 'Visible','off');
        set(handles.stepNum_label, 'Visible','off');
        set(handles.eps_edit, 'Visible','off');
        set(handles.eps_label, 'Visible','off');
        set(handles.workPrec_edit, 'Visible','off');
        set(handles.workPrec_label, 'Visible','off');
end


% --- Executes on selection change in methodGroupPopUp.
function methodGroupPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to methodGroupPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methodGroupPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methodGroupPopUp
methodGroup = handles.GuiHelpers.getCurrentPopupString(handles.methodGroupPopUp);
%set default methods for each group

%if strcmp(methodGroup,'LineSearch') == 0
    methods=strcat('Methods', filesep, 'MultiDimensional', filesep, methodGroup);
    multMethods=dir(methods);
    set(handles.multiDimMethodPopUp, 'Value', 1);
    set(handles.multiDimMethodPopUp, 'String', handles.GuiHelpers.fNames(multMethods));
    defaultMethodPos = handles.GuiHelpers.getDefaultMethodPosition(handles, methodGroup);
    set(handles.multiDimMethodPopUp, 'Value', defaultMethodPos);
    %multiDimMethodPopUp_Callback(hObject, eventdata, handles);
    defaultLineSearchPos = handles.GuiHelpers.getDefaultLineSearchPos(handles, '', methodGroup);
    set(handles.lineSearchPopUp, 'Value', defaultLineSearchPos);
    multiDimMethodPopUp_Callback(hObject, eventdata, handles);
    method = handles.GuiHelpers.getCurrentPopupString(handles.methodGroupPopUp);
    
    %lineSearchMethod = handles.GuiHelpers.getCurrentPopupString(handles.lineSearchPopUp);
    %handles.GuiHelpers.setLineSearchParams(handles, lineSearchMethod);
    
    if strcmp(handles.GuiHelpers.enableLineSearch(handles, methodGroup), 'None') == 1 || ...
       strcmp(handles.GuiHelpers.enableLineSearch(handles, method), 'None') == 1
        %set(handles.lineSearchPopUp, 'Visible', 'Off');
        %set(handles.lineSearch_label, 'Visible', 'Off');
    else
        %set(handles.lineSearchPopUp, 'Visible', 'On');
        %set(handles.lineSearch_label, 'Visible', 'On');
    end
    
    if ~handles.GuiHelpers.enableAdvancedPanel(handles, methodGroup) || ...
       ~handles.GuiHelpers.enableAdvancedPanel(handles, method)
        set(handles.cetiriPromenljive_panel, 'Visible', 'Off');
    else
        set(handles.cetiriPromenljive_panel, 'Visible', 'On');
    end
    
    set(handles.stop_cond_panel, 'Visible', 'On');
    
%end

% --- Executes during object creation, after setting all properties.
function methodGroupPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodGroupPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in oneDimButton.
function oneDimButton_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.oneDimPanel, 'Visible', 'On', 'Position',[5,8,1746,821]);
set(handles.twoDimPanel, 'Visible', 'Off');
set(handles.multDimPanel, 'Visible', 'Off');
% Hint: get(hObject,'Value') returns toggle state of oneDimButton


% --- Executes on button press in twoDim.
function twoDim_Callback(hObject, eventdata, handles)
% hObject    handle to twoDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.oneDimPanel, 'Visible', 'Off');
set(handles.twoDimPanel, 'Visible', 'On');
set(handles.multDimPanel, 'Visible', 'Off');
% Hint: get(hObject,'Value') returns toggle state of twoDim


% --- Executes on button press in multiDim.
function multiDim_Callback(hObject, eventdata, handles)
% hObject    handle to multiDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.oneDimPanel, 'Visible', 'Off');
set(handles.twoDimPanel, 'Visible', 'Off');
set(handles.multDimPanel, 'Visible', 'On');
% Hint: get(hObject,'Value') returns toggle state of multiDim


% --- Executes on selection change in oneDimMethodPopUp.
function oneDimMethodPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimMethodPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns oneDimMethodPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from oneDimMethodPopUp


% --- Executes during object creation, after setting all properties.
function oneDimMethodPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimMethodPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oneDimStartingPoint_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimStartingPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimStartingPoint as text
%        str2double(get(hObject,'String')) returns contents of oneDimStartingPoint as a double


% --- Executes during object creation, after setting all properties.
function oneDimStartingPoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimStartingPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oneDimStepSize_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimStepSize as text
%        str2double(get(hObject,'String')) returns contents of oneDimStepSize as a double


% --- Executes during object creation, after setting all properties.
function oneDimStepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oneDimMaxIter_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimMaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimMaxIter as text
%        str2double(get(hObject,'String')) returns contents of oneDimMaxIter as a double


% --- Executes during object creation, after setting all properties.
function oneDimMaxIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimMaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oneDimEps_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimEps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimEps as text
%        str2double(get(hObject,'String')) returns contents of oneDimEps as a double


% --- Executes during object creation, after setting all properties.
function oneDimEps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimEps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in oneDimFunctionPopUp.
function oneDimFunctionPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimFunctionPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns oneDimFunctionPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from oneDimFunctionPopUp


% --- Executes during object creation, after setting all properties.
function oneDimFunctionPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimFunctionPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oneDimGradEvalNum_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimGradEvalNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of oneDimGradEvalNum as text
%        str2double(get(hObject,'String')) returns contents of oneDimGradEvalNum as a double


% --- Executes during object creation, after setting all properties.
function oneDimGradEvalNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimGradEvalNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oneDimFEvalNum_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimFEvalNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimFEvalNum as text
%        str2double(get(hObject,'String')) returns contents of oneDimFEvalNum as a double


% --- Executes during object creation, after setting all properties.
function oneDimFEvalNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimFEvalNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oneDimCpuTime_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimCpuTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimCpuTime as text
%        str2double(get(hObject,'String')) returns contents of oneDimCpuTime as a double


% --- Executes during object creation, after setting all properties.
function oneDimCpuTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimCpuTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oneDimIterNum_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimIterNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimIterNum as text
%        str2double(get(hObject,'String')) returns contents of oneDimIterNum as a double


% --- Executes during object creation, after setting all properties.
function oneDimIterNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimIterNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oneDimFmin_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimFmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oneDimFmin as text
%        str2double(get(hObject,'String')) returns contents of oneDimFmin as a double


% --- Executes during object creation, after setting all properties.
function oneDimFmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oneDimFmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in oneDimMaximize.
function oneDimMaximize_Callback(hObject, eventdata, handles)
% hObject    handle to oneDimMaximize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%clear previous results
handles.GuiHelpers.clearOneDimResults(handles);
%load function and method
SELECTED_FUNCTION = handles.GuiHelpers.selectedFunction(handles.oneDimFunctionPopUp);
SELECTED_METHOD = handles.GuiHelpers.selectedMethod(handles.oneDimMethodPopUp);
%load method params
methodParams = loadOneDimMethodParams(handles);
%calculating
[ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = SELECTED_METHOD(SELECTED_FUNCTION, methodParams);
results = Results(fmin, xmin, valuesPerIter.gradientPerIteration(end), iterNum, cpuTime, evalNumbers, valuesPerIter);
%display results
handles.GuiHelpers.displayOneDimResults(results, handles);


% --- Executes on selection change in lineSearchPopUp.
function lineSearchPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to lineSearchPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lineSearchPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lineSearchPopUp
selectedLineSearchMethod = handles.GuiHelpers.getCurrentPopupString(handles.lineSearchPopUp);

if strcmp(handles.enabledLineSearch, 'All') || strcmp(selectedLineSearchMethod, handles.enabledLineSearch) == 1
    handles.lineSearchMethod = selectedLineSearchMethod;
    guidata(hObject, handles);
else
    currentPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), handles.lineSearchMethod));
    set(handles.lineSearchPopUp, 'Value', currentPos);
end
handles.GuiHelpers.setLineSearchParams(handles, selectedLineSearchMethod);
%starting step size needed only in this case(s)


% --- Executes during object creation, after setting all properties.
function lineSearchPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lineSearchPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function hesEvalNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to hesEvalNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hesEvalNumEdit as text
%        str2double(get(hObject,'String')) returns contents of hesEvalNumEdit as a double


% --- Executes during object creation, after setting all properties.
function hesEvalNumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hesEvalNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xminResultEdit_Callback(hObject, eventdata, handles)
% hObject    handle to xminResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xminResultEdit as text
%        str2double(get(hObject,'String')) returns contents of xminResultEdit as a double


% --- Executes during object creation, after setting all properties.
function xminResultEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xminResultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function outputPanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in defaultModeCheckbox.
function defaultModeCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to defaultModeCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of defaultModeCheckbox



function workPrec_edit_Callback(hObject, eventdata, handles)
% hObject    handle to workPrec_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of workPrec_edit as text
%        str2double(get(hObject,'String')) returns contents of workPrec_edit as a double


% --- Executes during object creation, after setting all properties.
function workPrec_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to workPrec_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ksi_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ksi_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ksi_edit as text
%        str2double(get(hObject,'String')) returns contents of ksi_edit as a double


% --- Executes during object creation, after setting all properties.
function ksi_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ksi_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function M_edit_Callback(hObject, eventdata, handles)
% hObject    handle to M_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_edit as text
%        str2double(get(hObject,'String')) returns contents of M_edit as a double


% --- Executes during object creation, after setting all properties.
function M_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
