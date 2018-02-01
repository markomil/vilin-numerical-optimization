function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = NewMethodTemplate( functionName, methodParams )
% ================================================================================================================
% Template for adding new optimization methods to Vilin. This is some most common example, method can be 
% implemented in any way as long as it returns correct values. :) To add new method change code in this file and 
% save in appropriate location, e.g. 'Methods/MultiDimensional/<some group>'.

% Code for reading parameter values is given at the begining.
% After that is shown how to evaluate function and set evalNumbers and valuesPerIter.
% In main loop is given code for creating params for line search, and evaluating selected line search method.
% Also at the end of each iteration valuesPerIter is updated.
% ================================================================================================================

    % fancy GUI message, remove when implementng method
    throw (MException('NumOpt:implementationError', 'Method %s not implemented.', mfilename));

    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    workPrec = methodParams.workPrec;
    xmin = methodParams.starting_point;
    t = methodParams.startingPoint;
    tPrev = t;
    it = 1;
    
    tic; % for measuring execution time
    % compute values for first iteration
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    % function value and gradient are computed, evalNumbers must be updated
    evalNumbers.incrementBy([1 1 0]); 
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    
    fPrev = fCurr + 1;

    % main method loop
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk =  % ... Determine direction vector
        % take vector of function values after each iteration
        fValues = valuesPerIter.functionPerIteration(1:it); 
        
        % creating params for line search method, 'dk' is search direction, and should be precomputed
        params = LineSearchParams(methodParams, fValues, grad, dk', xmin, tPrev, it);
        
        % executing line search method
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        
        % updating evaluation numbers, each line search can compute function value/graient/hessian multiple times
        evalNumbers = evalNumbers + lineSearchEvalNumbers; 
        
        % compute funcion and gradient value in current point xmin
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
	
        % method update code ... 
       
        % increment number of iterations and update values for last iteration
        it = it + 1;
        tPrev = t; % last value of step size t
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    fmin = fCurr;
    % getting elapsed time
    cpuTime = toc;
    % using only values for 'it' iterations
    valuesPerIter.trim(it);
    it = it - 1;
end

