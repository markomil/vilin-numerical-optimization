function [ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = RMSProp( functionName, methodParams )
% This is RMSProp (Root mean square Propagation) method with ineaxact line search rule

    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = methodParams.starting_point;
    t = methodParams.startingPoint;
            
    iterNum = 1;
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    valuesPerIter.setFunctionVal(iterNum, fCurr);
    valuesPerIter.setGradientVal(iterNum, norm(grad));
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    %numerical stability parameter
    ksi = 10^(-8); 
    % average decay of square gradients 
    adGrad = grad.^2;
    % root mean square of gradients
    rmsGrad = sqrt(adGrad + ksi);

    while (iterNum < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk = -1./rmsGrad.*grad;
        fValues = valuesPerIter.functionPerIteration(1:iterNum); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, dk', xmin, t, iterNum);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;

        fPrev = fCurr;        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        
        % update formula for average decay of square gradients vector 
        % according to the RMSProp update rule
        adGrad = 0.9*adGrad + 0.1*grad.^2;
        % root mean square of gradients
        rmsGrad = sqrt(adGrad + ksi);
        
        iterNum = iterNum + 1;
        valuesPerIter.setFunctionVal(iterNum, fCurr);
        valuesPerIter.setGradientVal(iterNum, norm(grad));
        valuesPerIter.setStepVal(iterNum, t);
    end

    fmin = fCurr;
    cpuTime = toc;
    valuesPerIter.trim(iterNum);
    iterNum = iterNum - 1;
end

