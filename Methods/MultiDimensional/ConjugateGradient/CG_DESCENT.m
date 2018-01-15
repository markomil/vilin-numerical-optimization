function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = CG_Descent( functionName, methodParams )
% ------------------      *******************          ------------------

% CG_Descent algorithm introduced by W.W. Hager and H. Zhang.
% see original paper: 
% W.W. Hager and H. Zhang. Algorithm 851: "CG_Descent, a conjugate 
% gradient method with guaranteed descent", 
% ACM Trans. Math. Software, 32(1):113-137, 2006.

% ------------------      *******************          ------------------


    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    starting_point = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;
	t = methodParams.startingPoint;
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
    D = 0.7;
    Q = 0;
    C = 0;
    w = 1e-3;
        
    methodParams.lineSearchMethod = @Wolfe;
    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, pk', xmin, t, it);
        
        Q = 1 + D * Q;
        C = C + ((fCurr - C) / Q);
        
        if (it >= 2 && abs(fValues(end) - fValues(end - 1)) <= w*C)
            methodParams.lineSearchMethod = @ApproxWolfe;
            params.ksi = 1e-6 * C;
        end
        
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
        betaCGD = (1 / py)*(yk - 2*pk*(norm(yk)^2 / py))'*grad;
        betaCGD = max(betaCGD, niK); % Restart
        pk = betaCGD*pk - grad;
        
        it = it + 1
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    cpuTime = toc;
    fmin = fCurr;
    valuesPerIter.trim(it);
    it = it - 1;
end
