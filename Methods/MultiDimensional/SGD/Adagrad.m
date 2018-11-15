function [ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = Adagrad( functionName, methodParams )
% This is adagrad (adaptive gradient) method with ineaxact line search rule

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
    % make vector of sum of square gradients
    sumGrad = grad.^2;
    % root of sum of square gradients
    rssGrad = sqrt(sumGrad + ksi);

    while (iterNum < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk = -1./rssGrad.*grad;
        fValues = valuesPerIter.functionPerIteration(1:iterNum); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, dk', xmin, t, iterNum);
        % update values
        fPrev = fCurr; 
        
        [t, xmin, fCurr, grad, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
       
        % update formula for sum of square gradients
        sumGrad = sumGrad + grad.^2;
        % root of sum of square gradients
        rssGrad = sqrt(sumGrad + ksi);
        
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

