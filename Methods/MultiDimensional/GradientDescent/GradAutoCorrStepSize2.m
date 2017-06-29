function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = GradAutoCorrStepSize2( funcionName, methodParams )
% Gradient descent method with step size auto correction 
% based on two previous iterations

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0, 0, 0);
    starting_point = methodParams.starting_point;
    delta = methodParams.step_size; 
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;
    d = delta;
    it = 1;

    [val, grad] = feval(funcionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    valuesPerIter.setFunctionVal(it, val);
    valuesPerIter.setGradientVal(it, norm(grad));
    valuesPerIter.setStepVal(it, d);
    it = it+1;
    fp2 = val;
    xmin = xmin - d * grad'/norm(grad);

    % compute val and grad after first iteration
    [val, grad] = feval(funcionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    valuesPerIter.setFunctionVal(it, val);
    valuesPerIter.setGradientVal(it, norm(grad));
    fp1 = val;
    xmin = xmin - d * grad'/norm(grad);

    % process
    while (it < maxIter && norm(grad) > epsilon)
        
        [val, grad, ~] = feval(funcionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);

        if val < fp1 && fp1 < fp2
            d = d * 1.5;
        else if val > fp1
                d = d * 0.5;
            end
        end
        fp2 = fp1;
        fp1 = val;
        xmin = xmin - d * grad'/norm(grad);

        it = it+1;
        valuesPerIter.setFunctionVal(it, val);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, d);
    end
    
    fmin = val;
    cpuTime = toc;
    valuesPerIter.trim(it);
    it = it - 1;
    
end
