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
    it = 1;
    
    tic; % for measuring execution time
    % compute values for first iteration
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]); % function value and hessian are computed, 
    				      % evalNumbers must be updated
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    
    fPrev = fCurr + 1;

    % main method loop
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)

        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it, grad, dk, methodParams.startingPoint); % 'lsStartPnt' is starting point for line search
        params = LineSearchParams(methodParams, fCurr, grad, dk', xmin, lsStartPnt); % creating params for line search method, 
										     % 'dk' is search direction, and should be precomputed
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params); % executing line search method
        evalNumbers = evalNumbers + lineSearchEvalNumbers; % updating evaluation numbers, 
							   % each line search can compute function value/graient/hessian multiple times

	% method update code ... 

 	% increment number of iterations and update values for last iteration
        it = it + 1;
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    fmin = fCurr;
    cpuTime = toc; % getting elapsed time
    valuesPerIter.trim(it); % using only values for 'it' iterations
    it = it - 1;
end

