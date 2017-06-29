function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = Levenberg( functionName, methodParams )
%%%%%%%%                Header              %%%%%%%%%%
%       This is Modified Newton method implemented by 
%       using numerical gradient and Hessian computations.
%       It is invented by Levenberg
%       No line search methods are needed for computing step 
%       size in every iteration.
%
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    xmin = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    doRecalculate = true;
    lambda = 10^(-4);
    lambdaMin = 10^(-10);
    lambdaMax = 10^(7);
    lamMul = 10;
    maxLambdaNotAchieved = true;
    dim = length(xmin);
    
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    [val, gr, ~] = feval(functionName,xmin,[1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, val);
    valuesPerIter.setGradientVal(it, grNorm);
                
    % process
    while  (it < maxIter)  && maxLambdaNotAchieved && (grNorm > eps)
        
        % Computes gradient and Hessian
        if doRecalculate
            if it > 1
                [~ , gr, Hes] = feval(functionName, xmin, [0 1 1]);   
                evalNumbers.incrementBy([0 1 1]);
            else
                [~ , ~, Hes] = feval(functionName, xmin, [0 1 1]);   
                evalNumbers.incrementBy([0 0 1]);
            end;
            grNorm = double(norm(gr));
        end
        
        % Computes search direction according to the Levenberg rule 
        leftTerm = Hes + lambda*diag(ones(dim,1));
        dir = -(leftTerm\gr)';
        
        params = LineSearchParams(methodParams, val, gr, dir, xmin, 1);
        [t, xminCurr, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        it = it + 1;
        
        % compute value in new point
        [valCurr , ~, ~] = feval(functionName,xminCurr, [1 0 0]);   
        evalNumbers.incrementBy([1 0 0]);
            
        if valCurr >= val
            doRecalculate = false;
            if lambda == lambdaMax
                maxLambdaNotAchieved = false;
	        else
                lambda = min(lamMul*lambda, lambdaMax);
            end
        else
            % update paramaters
            xmin = xminCurr;
            lambda = max(lambda/lamMul, lambdaMin);
            val = valCurr;
            doRecalculate = true;
        end;
                                  
        valuesPerIter.setFunctionVal(it, val);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
    end;
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    fmin = val;
    it = it - 1;
end
