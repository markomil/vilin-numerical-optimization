function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = GoldsteinPrice( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *      Goldstein Price method       *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Goldstein Price is Modified Newton method implemented by using 
%   numerical gradient and Hessian computations. The general idea is to 
%   use the steepest descent direction in situations when Hessian is not 
%   positive definite, and thus not guarantee the convergence. Parametar 
%   eta is a threshold, in so called angle rule, which determines wheater 
%   to use gradient or Newton search direction. Finally, appropriate line 
%   search is used for computing step size in every iteration.

%   A.A. Goldstein, J.F. Price,
%   An effective algorithm for minimization, 
%   Numer. Math. 10 (1967) 184–189.

%   ------------------      *******************        ------------------
  
    % set initial values
    tic;                                    % to compute CPU time
    evalNumbers = EvaluationNumbers(0,0,0);
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    eta = 0.2;                              % parametar eta initialization
    xmin = methodParams.starting_point;
    t = methodParams.startingPoint;
    it = 1;                                 % number of iteration
    
    [fCurr, grad, Hes] = feval(functionName, xmin, [1 1 1]);
    evalNumbers.incrementBy([1 1 1]);
    grNorm = double(norm(grad));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    % add values for plot
    if (size(xmin, 2) == 2)
        valuesPerIter.setXVal(it, xmin);
    end
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
                
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        % Computes dir according to the Goldstein Price rule 
        dir = -Hes\grad;                  % computes Newton direction
        if dir'*(-grad)/(norm(dir)*norm(grad)) < eta || sum(isnan(dir)) > 0
            dir = -grad;
        end
                
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, dir', xmin, t, it);
        % update values
        fPrev = fCurr;
                
        % Computes xmin and step-size according to the line search method rule
        [t, xmin, fCurr, grad, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        grNorm = double(norm(grad));
                
        % compute numerical Hessian in new point
        [~ , ~, Hes] = feval(functionName, xmin, [0 0 1]);
        evalNumbers.incrementBy([0 0 1]);
                
        it = it + 1;
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
        % add values for plot
        if (size(xmin, 2) == 2)
            valuesPerIter.setXVal(it, xmin);
            valuesPerIter.setDirVal(it, dir);
        end
    end
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    fmin = fCurr;
    it = it - 1;
end
