function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = CG_DESCENT( functionName, methodParams )
    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    starting_point = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;
    it = 1;
    ni = 0.01;
    
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
        
        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it+1, grad, pk, methodParams.startingPoint);
        params = LineSearchParams(methodParams, fCurr, grad, pk', xmin, lsStartPnt);
        [t, xmin, lineSearchEvalNumbers] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        % update values
        fPrev = fCurr;
        gradOld = grad;
        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);

        % compute parameter beta
        niK = -1 / (norm(pk) * min(ni, norm(grad)));
        yk = grad - gradOld;
        py = pk' * yk;
        betaCGD = (1 / py)*(yk - 2*pk*(norm(yk)^2 / py))'*grad; %  (grad'*(grad-gradOld))/(gradOld'*gradOld);
        betaCGD = max(betaCGD, niK); % Restart
        pk = betaCGD*pk - grad;
        
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
