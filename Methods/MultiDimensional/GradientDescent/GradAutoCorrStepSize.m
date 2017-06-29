function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = GradAutoCorrStepSize( functionName, methodParams )
% Gradient descent method with step size auto correction 
% based on previous iteration

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    starting_point = methodParams.starting_point;
    delta = methodParams.step_size; 
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;
    d = delta;
    it = 1;

    [val, grad] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    valuesPerIter.setFunctionVal(it, val);
    valuesPerIter.setGradientVal(it, norm(grad));
    valuesPerIter.setStepVal(it, d);
    fp1 = val;
    xmin = xmin - d * grad'/norm(grad);

    % process
    while (it < maxIter && norm(grad) > epsilon)
        
        [val,grad] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        if val > fp1
            d = d * 0.5;
        end
        fp1 = val;
        xmin = xmin - d * grad'/norm(grad);

        it = it + 1;
        valuesPerIter.setFunctionVal(it, val);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, d);
    end
    
    fmin = val;
    cpuTime = toc;
    valuesPerIter.trim(it);
    it = it - 1;
end
