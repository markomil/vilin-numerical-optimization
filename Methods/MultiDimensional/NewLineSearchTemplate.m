function [ outT, outX, evalNumbers ] = NewLineSearchTemplate( functionName, params )
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
    val = params.val;
    gr = params.grad;
    dir = params.dir;
    rho = params.rho;
    it = 1; % number of iteration
    t1 = params.tStart; % starting value for t
    tInit = params.tPrev; % value of step t from previous iteration of the main method
    
    [val1,~] = feval(functionName, x0+t1*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
          
    % main loop  
    while 0 - %line search condition, e.g. Armijo, Goldstein, Wolfe...
        
        % method code
        % should contain evalNumbers.incrementBy(...)
        it = it + 1;
    end 
end
