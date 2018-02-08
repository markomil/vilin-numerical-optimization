function p = presets()

p = struct('getMethod', @getMethod,...
           'getLineSearch', @getLineSearch,...
           'getAllowedLineSearch', @getAllowedLineSearch,...
           'getLineSearchParams', @getLineSearchParams);
             
function method = getMethod(methodGroup)
    method = ''; % e.g. current
    switch methodGroup
        case 'ConjugateGradient'
            method = 'PolakRibiere';
        case 'GradientDescent'
            method = 'GradientLineSearch';
        otherwise
            
    end
    

function lineSearch = getLineSearch(method, methodGroup)
    lineSearch = ''; % e.g. current
    if strcmp(methodGroup,'ConjugateGradient') == 1
        lineSearch = 'StrongWolfe';
    end
    switch method
        case 'GoldsteinPrice'
            lineSearch = 'Wolfe';
        case 'Levenberg'
            lineSearch = 'FixedStepSize';
        case 'LevenbergMarquardt'
            lineSearch = 'FixedStepSize';
        case 'GradientLineSearch'
            lineSearch = 'Wolfe';
        case 'GradientAccelerated'
            lineSearch = 'Wolfe';
        case 'BarzilaiBorwein'
            lineSearch = 'NonMonotone';
        case 'ScalarCorrection'
            lineSearch = 'NonMonotone';
        case 'BFGS'
            lineSearch = 'Wolfe';
        case 'DFP'
            lineSearch = 'Wolfe';
        case 'NewtonLineSearch'
            lineSearch = 'FixedStepSize';
        case 'CG_Descent'
            lineSearch = 'ApproxWolfe';
        otherwise
            
    end
    
function allowedLineSearch = getAllowedLineSearch(method)
    switch method
        case 'Levenberg'
            allowedLineSearch = 'None';
        case 'LevenbergMarquardt'
            allowedLineSearch = 'None';
        case 'TrustRegion'
            allowedLineSearch = 'None';
        otherwise
            allowedLineSearch = 'All';
    end

function params = getLineSearchParams(lineSearch)
    beta = 0.8;
    startingPoint = 1;
    M = 10;
    sigma = 0.01;
    rho = 0.001;
    theta = 0.5; %  used in the update rule in ApproxWolfe
    gamma = 0.66; % determines when a bisection step is performed in ApproxWolfe
    w = 1e-3;
    
    if strcmp(lineSearch, 'ApproxWolfe') == 1
        sigma = 0.9;
        rho = 0.1;
    end
    
    params = struct('beta', beta, 'starting_point', startingPoint, 'M', M, ...
        'sigma', sigma, 'rho', rho, 'theta', theta, 'gamma', gamma, 'w', w);