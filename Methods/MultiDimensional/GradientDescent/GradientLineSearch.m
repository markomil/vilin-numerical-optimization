function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = GradientLineSearch( functionName, methodParams )
% Gradient descent method with ineaxact line search rule

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = methodParams.starting_point;
    it = 1;
    
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;

    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk = -grad;
        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it, grad, dk, methodParams.startingPoint);
        params = LineSearchParams(methodParams, fCurr, grad, dk', xmin, lsStartPnt);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;

        fPrev = fCurr;        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);

        it = it + 1;
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    fmin = fCurr;
    cpuTime = toc;
    valuesPerIter.trim(it);
    it = it - 1;
end

