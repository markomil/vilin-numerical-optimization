function [ fmin, xmin, iterNum, cpuTime, evalNumbers, valuesPerIter ] = Adadelta( functionName, methodParams )
% This is Adadelta method with ineaxact line search rule

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
    gamma = 0.9;
    
    n = length(xmin);
    %numerical stability parameter
    ksi = 10^(-8); 
    % average decay of square parameters update 
    adUpdate = zeros(n,1);
    % average decay of square gradients 
    adGrad = grad.^2;
    % root mean square of gradients
    rmsGrad = sqrt(adGrad + ksi);
    % root mean square of parameters update
    rmsUpdate = ones(n,1);
    
    while (iterNum < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk = -(rmsUpdate./rmsGrad).*grad;
        fValues = valuesPerIter.functionPerIteration(1:iterNum); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, dk', xmin, t, iterNum);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;

        fPrev = fCurr;        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        
        % update formula for average decay of square gradients 
        adGrad = gamma*adGrad + (1-gamma)*grad.^2;  
        % root mean square of gradients
        rmsGrad = sqrt(adGrad + ksi);
        % current square parameters update 
        currUpdate = (t*dk).^2;
        % update formula for average decay of square parameters update 
        adUpdate = gamma*adUpdate + (1-gamma)*currUpdate;
        % root mean square of parameters update
        rmsUpdate = sqrt(adUpdate + ksi);
            
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

