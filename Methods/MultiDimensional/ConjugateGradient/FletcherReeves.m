function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = FletcherReeves( functionName, methodParams )
% Fletcher-Reeves version of Conjugate gradient method

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    starting_point = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;

    it = 1;
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    pk = - grad;
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it, grad, pk, methodParams.startingPoint);
        params = LineSearchParams(methodParams, fCurr, grad, pk', xmin, lsStartPnt);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        % update values
        fPrev = fCurr;
        gradOld = grad;
        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        
        % compute parameter beta
        betaFR = (grad'*grad)/(gradOld'*gradOld);
        %betaFR = max(betaFR, 0); % Restart
        pk = betaFR*pk - grad;
        
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
