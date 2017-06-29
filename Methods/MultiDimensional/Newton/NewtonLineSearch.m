function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = NewtonLineSearch( functionName, methodParams )
% Ordinary Newton's method with appropriate line search rule

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0, 0, 0);
    starting_point = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;
    it = 1;
    
    [fCurr, grad, hes] = feval(functionName, xmin, [1 1 1]);
    evalNumbers.incrementBy([1 1 1]);
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;

    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        % Computes Newton search direction
        dk = -hes\grad;
        
        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it, grad, dk, methodParams.startingPoint);
        params = LineSearchParams(methodParams, fCurr, grad, dk', xmin, lsStartPnt);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
              
        fPrev = fCurr;
        [fCurr, grad, hes] = feval(functionName, xmin, [1 1 1]);
        evalNumbers.incrementBy([1 1 1]);

        it = it + 1;
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    cpuTime = toc;
    fmin = fCurr;
    valuesPerIter.trim(it);
    it = it - 1;
end

