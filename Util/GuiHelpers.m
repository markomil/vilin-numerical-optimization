function functions = GuiHelpers()
functions = struct('displayMultDimResults',@displayMultDimResults, ...
                  'displayOneDimResults', @displayOneDimResults,...
                  'clearMultDimResults', @clearMultDimResults,...
                  'clearOneDimResults', @clearOneDimResults,...
                  'plotMultDimResults', @plotMultDimResults,...
                  'loadMultDimMethodParams', @loadMultDimMethodParams,...
                  'loadOneDimMethodParams', @loadOneDimMethodParams,...
                  'updateStartingPoint', @updateStartingPoint,...
                  'startingPointValid', @startingPointValid,...
                  'initSliders', @initSliders,...
                  'fNames', @fNames,...
                  'subDirs', @subDirs,...
                  'setGraphicTitles', @setGraphicTitles,...
                  'getCurrentPopupString', @getCurrentPopupString,...
                  'prepParams', @prepParams,...
                  'getVector', @getVector,...
                  'selectedFunction', @selectedFunction,...
                  'selectedMethod', @selectedMethod,...
                  'startingPoint', @startingPoint,...
                  'setStartingPoint', @setStartingPoint,...
                  'stepSize', @stepSize,...
                  'dimension', @dimension,...
                  'setDimension', @setDimension,...
                  'numberOfSteps', @numberOfSteps,...
                  'eps', @eps,...
                  'delta_min', @delta_min,...
                  'od', @od,...
                  'do', @do,...
                  'beta', @beta,...
                  'sigma', @sigma,...
                  'rho', @rho,...
                  'startingPointZaCetiriPromenljive', @startingPointZaCetiriPromenljive,...
                  'clearGraphics', @clearGraphics,...
                  'plotGrad', @plotGrad,...
                  'plotVal', @plotVal,...
                  'crtaj_panel', @crtaj_panel,...
                  'pressCtrl_C', @pressCtrl_C,...
                  'releaseCtrl_C', @releaseCtrl_C,...
                  'getDefaultMethodPosition', @getDefaultMethodPosition,...
                  'getDefaultLineSearchPos', @getDefaultLineSearchPos,...
                  'isDefaultMode', @isDefaultMode,...
                  'enableLineSearch', @enableLineSearch);

function displayMultDimResults(results, handles)
set(handles.resultMethodName, 'String', getCurrentPopupString(handles.multiDimMethodPopUp));
set(handles.resultFunctionName, 'String', getCurrentPopupString(handles.multFunctionPopUp));
set(handles.fminResultEdit, 'String',results.fmin);
len = 100; 
if length(results.xmin) < 100
    len = length(results.xmin);
end;
set(handles.xminResultEdit, 'String',strrep(strrep(mat2str(results.xmin(1:len), 3), '[', ''),']',''));
set(handles.iterNumResultEdit, 'String', results.iterNum);
set(handles.gradNormResultEdit, 'String', results.gradNorm);
set(handles.cpuTimeEdit, 'String', results.cpuTime);
set(handles.funEvalNumEdit, 'String', results.evalNumbers.functionEvalNo);
set(handles.gradEvalNumEdit, 'String', results.evalNumbers.gradientEvalNo);
set(handles.hesEvalNumEdit, 'String', results.evalNumbers.hessianEvalNo);

function displayOneDimResults(results, handles)
set(handles.oneDimFmin, 'String', results.fmin);
set(handles.oneDimIterNum, 'String', results.iterNum);
set(handles.oneDimCpuTime, 'String', results.cpuTime);
set(handles.oneDimFEvalNum, 'String', results.evalNumbers.functionEvalNo);

function clearMultDimResults(handles)
set(handles.resultMethodName, 'String', '');
set(handles.resultFunctionName, 'String', '');
set(handles.fminResultEdit, 'String','');
set(handles.xminResultEdit, 'String','');
set(handles.iterNumResultEdit, 'String', '');
set(handles.cpuTimeEdit, 'String', '');
set(handles.funEvalNumEdit, 'String', '');
set(handles.gradEvalNumEdit, 'String', '');
set(handles.hesEvalNumEdit, 'String', '');
set(handles.gradNormResultEdit, 'String', '');

function clearOneDimResults(handles)
set(handles.oneDimFmin, 'String', '');
set(handles.oneDimIterNum, 'String', '');
set(handles.oneDimCpuTime, 'String', '');
set(handles.oneDimFEvalNum, 'String', '');

function plotMultDimResults(results, handles)
plot(handles.axesGr, results.valuesPerIteration.iterations, results.valuesPerIteration.gradientPerIteration);
plot(handles.axesVal, results.valuesPerIteration.iterations, results.valuesPerIteration.functionPerIteration);
setGraphicTitles(handles);

function [methodParams, success, message] = loadMultDimMethodParams(handles,functionName)
success = true;
message = '';
STARTING_POINT = startingPoint(handles.startingPointEdit);
STEP_SIZE = stepSize(handles.stepSizeEdit);
DIMENSION = dimension(handles.dimensionEdit);
NUMBER_OF_STEPS = numberOfSteps(handles.stepNumEdit);
DELTA_MIN = delta_min(handles.deltaMin_edit);
OD = od(handles.od_edit);
DO = do(handles.do_edit);
EPS = eps(handles.epsEdit);
WORKPREC = workPrecision(handles.workPrecEdit);
BETA = beta(handles.beta_edit);
SIGMA = sigma(handles.sigma_edit);
RHO = rho(handles.rho_edit);
KSI = ksi(handles.ksi_edit);
EM = M(handles.M_edit);
STARTINGPOINT = startingPoint(handles.startingPoint_edit);
if DIMENSION ~= length(STARTING_POINT)
    [STARTING_POINT, validSP, message] = updateStartingPoint(handles, functionName, DIMENSION);
    success = success && validSP;
end
LINE_SEARCH_METHOD = str2func(handles.lineSearchMethod);
methodParams  = MathodParams(STARTING_POINT, STEP_SIZE, DIMENSION, NUMBER_OF_STEPS, EPS, WORKPREC, DELTA_MIN, OD, DO, BETA, SIGMA, RHO, EM, KSI, STARTINGPOINT, LINE_SEARCH_METHOD);

function methodParams = loadOneDimMethodParams(handles)
STARTING_POINT = startingPoint(handles.oneDimStartingPoint);
STEP_SIZE = stepSize(handles.oneDimStepSize);
DIMENSION = 1;
NUMBER_OF_STEPS = numberOfSteps(handles.oneDimMaxIter);
DELTA_MIN = delta_min(handles.oneDimEps);
OD = 0;
DO = 1;
EPS = eps(handles.oneDimEps);
BETA = 1;
SIGMA = 1;
RHO = 1;
STARTINGPOINT_zaCetiriPromenljive = startingPointZaCetiriPromenljive(handles.startingPoint_edit);
methodParams  = MathodParams(STARTING_POINT, STEP_SIZE, DIMENSION, NUMBER_OF_STEPS, EPS, DELTA_MIN, OD, DO, BETA, SIGMA, RHO, STARTINGPOINT_zaCetiriPromenljive, 'null');

function [newSP, valid, message] = updateStartingPoint(handles, functionName, newDimension)
%returns new starting point and updates it's points value in GUI
sps = StartingPointGenerator(newDimension);
newSP = sps(functionName);
clear sps;
[valid, message] = startingPointValid(functionName, newSP);
if valid
    handles.startingPoints(functionName) = newSP;
    setStartingPoint(handles, newSP);
else
    newSP = [];
end

function [valid, message] = startingPointValid(functionName, startingPoint)
switch functionName
    case 'GenRosenbrock'
        valid = length(startingPoint) > 1;
        message = 'Starting point dimension must be 2 or greater.';
    otherwise
        valid = true;
        message = 'If you see this there is a bug in code.';
end

function initSliders(handles, length)
set(handles.gradEndSlider, 'Max', length - 1);
set(handles.gradEndSlider, 'Min', min(2, length-1));
set(handles.gradEndSlider, 'Value', length - 1);
set(handles.gradStartSlider, 'Max', length - 1);
set(handles.gradStartSlider, 'Min', 0);
set(handles.gradStartSlider, 'Value', 0);
%set values in boxes
set(handles.gradPlotStart, 'Value', 0);
set(handles.gradPlotStart, 'String', 0);
set(handles.gradPlotEnd, 'Value', length - 1);
set(handles.gradPlotEnd, 'String', length - 1);

function namess = fNames(dir)
names = regexp({dir(3:end).name}, '(\w*.m$)','match');
names = [names{:}];
    for i=1:length(names)
        names(i) = strtok(names(i), '.');
    end
namess = reshape(names, [length(names),1]);

function subFolders=subDirs(mainDir)
files = mainDir(3:end);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = {files(dirFlags).name};

function setGraphicTitles(handles)
title(handles.axesGr, 'Gradient value per iteration');
xlabel(handles.axesGr, 'Iteration');
ylabel(handles.axesGr, 'Gradient value');
%-------------------------
title(handles.axesVal, 'Function value per iteration');
xlabel(handles.axesVal, 'Iteration');
ylabel(handles.axesVal, 'Function value');

function str = getCurrentPopupString(hh)
%# getCurrentPopupString returns the currently selected string in the popupmenu with handle hh
if ~ishandle(hh) || strcmp(get(hh,'Type'),'popupmenu')
error('getCurrentPopupString needs a handle to a popupmenu as input')
end

%# get the string - do it the readable way
list = get(hh,'String');
val = get(hh,'Value');
if iscell(list)
   str = list{val};
else
   str = list(val,:);
end

function prepParams(handles)
sf = selectedFunction(handles.multFunctionPopUp);
if handles.startingPoints.isKey(sf)
    setStartingPoint(handles, handles.startingPoints(sf));
    setDimension(handles, length(handles.startingPoints(sf)));
else
    setStartingPoint(handles, 1);
    setDimension(handles, 1);
end

function vector = getVector(hh)
% Gets vector from edit textbox
str=get(hh,'String');
tmp = regexp (str,',','split');
vector=cellfun(@str2double,tmp);

function sp = startingPoint(startingPointEdit)
sp = getVector(startingPointEdit);

function sf = selectedFunction(functionPopUP)
sf = strtok(getCurrentPopupString(functionPopUP),'.');

function sm = selectedMethod(methodPopUp)
sm = str2func(strtok(getCurrentPopupString(methodPopUp),'.'));

function setStartingPoint(handles, sp)
str = strrep(strrep(strrep(mat2str(sp),']',''),'[',''),' ',',');
set(handles.startingPointEdit,'String',str);

function ss = stepSize(stepSizeEdit)
ss=getVector(stepSizeEdit);

function d = dimension(dimensionEdit)
d = getVector(dimensionEdit);

function setDimension(handles, d)
set(handles.dimensionEdit,'String',d);

function nos = numberOfSteps(stepNumEdit)
nos = getVector(stepNumEdit);

function e = eps(epsEdit)
e = getVector(epsEdit);

function wp = workPrecision(wpEdit)
    wp = getVector(wpEdit);

function delta_min = delta_min(deltaMin_edit)
delta_min = getVector(deltaMin_edit);

function od = od(od_edit)
od = getVector(od_edit);

function do = do(do_edit)
do = getVector(do_edit);

function beta = beta(beta_edit)
beta = getVector(beta_edit);

function sigma = sigma(sigma_edit)
sigma = getVector(sigma_edit);

function rho = rho(rho_edit)
rho = getVector(rho_edit);

function ksi = ksi(ksi_edit)
ksi = getVector(ksi_edit);

function M = M(M_edit)
M = getVector(M_edit);

function startingPoint = startingPointZaCetiriPromenljive(startingPoint_edit)
startingPoint = getVector(startingPoint_edit);

function clearGraphics(handles)
plot(handles.axesGr, [0], [0]);
plot(handles.axesVal, [0], [0]);
setGraphicTitles(handles);

function plotGrad(handles)
plotStart = getVector(handles.gradPlotStart) + 1;
plotEnd = getVector(handles.gradPlotEnd) + 1;
plot(handles.axesGr, handles.iterations(plotStart:plotEnd), handles.gradPerIter(plotStart:plotEnd));
plot(handles.axesVal, handles.iterations(plotStart:plotEnd), handles.valuesPerIter(plotStart:plotEnd));
setGraphicTitles(handles);

function plotVal(handles)
plotStart = getVector(handles.valPlotStart) + 1;
plotEnd = getVector(handles.valPlotEnd) + 1;
plot(handles.axesVal, handles.iterations(plotStart:plotEnd), handles.valuesPerIter(plotStart:plotEnd));
setGraphicTitles(handles);

%nikola
function crtaj_panel(handles)
set(handles.deltaMin_panel, 'Visible','off');
set(handles.constKorak_panel, 'Visible','off');

function pressCtrl_C
    import java.awt.Robot;
    import java.awt.event.*;
    SimKey=Robot;
    SimKey.keyPress(KeyEvent.VK_CONTROL);
    SimKey.keyPress(KeyEvent.VK_C);
    
function releaseCtrl_C(ignore1, ignore2)
    import java.awt.Robot;
    import java.awt.event.*;
    SimKey=Robot;
    SimKey.keyRelease(KeyEvent.VK_CONTROL);
    SimKey.keyRelease(KeyEvent.VK_C);
    
function defaultMethodPos = getDefaultMethodPosition(handles, methodGroup)
    switch methodGroup
        case 'ConjugateGradient'
            defaultMethodPos = find(strcmp(get(handles.multiDimMethodPopUp, 'String'), 'PolakRibiere'));
        case 'GradientDescent'
            defaultMethodPos = find(strcmp(get(handles.multiDimMethodPopUp, 'String'), 'GradientLineSearch'));
        otherwise
            defaultMethodPos = 1;
    end
    
function defaultLineSearchPos = getDefaultLineSearchPos(handles, method, methodGroup)
if ~ isDefaultMode(handles)
    defaultLineSearchPos = get(handles.lineSearchPopUp, 'Value');
else
    switch methodGroup
        case 'ConjugateGradient'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'StrongWolfe'));
        otherwise
            defaultLineSearchPos = get(handles.lineSearchPopUp, 'Value');
    end
    switch method
        case 'GoldsteinPrice'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'Wolfe'));
        case 'Levenberg'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'FixedStepSize'));
        case 'LevenbergMarquardt'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'FixedStepSize'));
        case 'GradientLineSearch'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'Wolfe'));
        case 'GradientAccelerated'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'Wolfe'));
        case 'BarzilaiBorwein'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'NonMonotone'));
        case 'ScalarCorrection'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'NonMonotone'));
        case 'BFGS'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'Wolfe'));
        case 'DFP'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'Wolfe'));
        case 'NewtonLineSearch'
            defaultLineSearchPos = find(strcmp(get(handles.lineSearchPopUp, 'String'), 'FixedStepSize'));
        otherwise
            
    end
end

    
function isDefault = isDefaultMode(handles)
    isDefault = (get(handles.defaultModeCheckbox, 'Value') == 1.0);


    function enabledLineSearch = enableLineSearch(handles, method)
    switch method
        case 'Levenberg'
            enabledLineSearch = 'FixedStepSize';
        case 'LevenbergMarquardt'
            enabledLineSearch = 'FixedStepSize';
        otherwise
            enabledLineSearch = 'All';
    end
