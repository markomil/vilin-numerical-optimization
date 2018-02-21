function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = GradientLineSearch( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *       Gradient Line Search        *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Gradient line search method is classical gradient descent (steepest
%   descent) method, originally developed by M. A. Cauchy for solving the 
%   system of linear equation. One of the most famous and popular method
%   which can be used for solving unconstrained optimization nonlinear 
%   problems after appliying appropriate line search procedure. 

%   M. A. Cauchy,
%   Methode generale pour la resolution des systemes d'equations simultanees, 
%   Comp. Rend. Acad. Sci. Par., 25 (1847), 536--538.

%   ------------------      *******************        ------------------

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0, 0, 0);
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = methodParams.starting_point;
    t = methodParams.startingPoint;
    it = 1;
    
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;

    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dk = -grad;
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, dk', xmin, t, it);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;

        fPrev = fCurr;        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);

        it = it + 1;
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    fmin = fCurr;
    cpuTime = toc;
    valuesPerIter.trim(it);
    it = it - 1;
end

