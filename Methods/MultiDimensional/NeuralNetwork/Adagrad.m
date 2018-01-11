function [ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = Adagrad( functionName, methodParams )
% This is adagrad (adaptive gradient) method with ineaxact line search rule
% Default line search is fixed line search with appropriatelly chosen value. 
% Let's say that the line search (learning rate) value depend on the nature
% of the goal function that we are minimizing. We can put default lr = 0.01;

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
    % make vector of lerning rates for each unknown parameter
    G = grad.^2;
    learRate = 1./sqrt(G + ksi);

    while (iterNum < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk = -learRate.*grad;
        fValues = valuesPerIter.functionPerIteration(1:iterNum); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, dk', xmin, t, iterNum);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;

        fPrev = fCurr;        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        
        % update the vector of lerning rates 
        % according to the adagrad update rule
        G = G + grad.^2;
        learRate = 1./sqrt(G + ksi);
        
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

