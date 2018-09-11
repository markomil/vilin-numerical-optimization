function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = Levenberg( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *         Levenberg  method         *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Levenberg  method is s originally constructed by Levenberg and is 
%   well studied as a solver for non-linear least squares problems. 
%   In order to use it as a solver for unconstrained optimization problems
%   some small changes are done. Namely, instead of J^TJ which represent 
%   the Hessian approximation (for least squares problems) the true Hessian 
%   can be used in the case of minimizing general nonlinear functions.
%   Trust region strategy is imposed ad thus no line search methods are 
%   needed for computing step size in every iteration.

%   K. Levenberg,
%   Method for the Solution of Certain Non-Linear Problems in Least Squares, 
%   Quarterly of Applied Mathematics, 2 (1944) 164-168.

%   ------------------      *******************        ------------------

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
    t = methodParams.startingPoint;
    
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
            end
            grNorm = double(norm(gr));
        end
        
        % Computes search direction according to the Levenberg rule 
        leftTerm = Hes + lambda*diag(ones(dim,1));
        dir = -(leftTerm\gr)';
        
        % computes next point with fixed step size 
        % methodParams.lineSearchMethod = 'FixedStepSize';
        xminCurr = xmin + t*dir;
                        
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
        end
        
        it = it + 1;
                                  
        valuesPerIter.setFunctionVal(it, val);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
    end
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    fmin = val;
    it = it - 1;
end
