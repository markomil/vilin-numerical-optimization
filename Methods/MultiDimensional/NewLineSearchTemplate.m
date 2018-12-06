function [ outT, outX, outVal, outGr, evalNumbers ] = NewLineSearchTemplate( functionName, params )
% ================================================================================================================
% Template for adding new line search methods to Vilin. This is some most common example, method can be 
% implemented in any way as long as it returns correct values. :) To add new line search method change code in this 
% file and save in 'Methods/MultiDimensional/LineSearch.

% Code for reading parameter values is given at the begining.
% After that is shown how to evaluate function and set evalNumbers.
% Main loop should contain exit condition and method's code.
% ================================================================================================================

    % fancy GUI message, remove when implementng method
    throw (MException('NumOpt:implementationError', 'Line search method %s not implemented.', mfilename));
    
    % read initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val0 = vals(end); % take last (current) function value
    gr = params.grad;
    dir = params.dir;
    rho = params.rho;
    sigma = params.sigma;
    ksi = params.ksi;
    tInit = params.tInitStart; % starting value for t
    it = 1; % number of iteration
    tPrev = params.tPrev; % value of step t from previous iteration of the main method
    
    [val1, gr1, ~] = feval(functionName, x0+tInit*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
          
    % main loop  
    while 0 - %line search condition, e.g. Armijo, Goldstein, Wolfe...
        
        % method code
        % should contain evalNumbers.incrementBy(...)
        it = it + 1;
    end 
    
    % save output values
    outX = x0 + t*dir;
    outT = t;
    outVal = val1;
    outGr = gr1;
    
    % if gradient is not computed in current point x0 + t*dir then
    % use the following code
        [~, outGr, ~] = feval(functionName, outX, [0 1 0]);   
        evalNumbers.incrementBy([0 1 0]);
end
