function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = FletcherReeves( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *      Fletcher-Reeves method       *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   This is nonlinear conjugate gradient method for solving 
%   large-scale unconstrained minimization problem originally designed 
%   by Fletcher and Reeves. In order to converge Strong Wolfe line 
%   search should be applied.

%   R. Fletcher, C.M. Reeves,
%   Function minimization by conjugate gradients, 
%   Comput. J. 7 (1964) 149–154.

%   ------------------      *******************        ------------------

    % set initial values
    tic;
    evalNumbers = EvaluationNumbers(0,0,0);
    starting_point = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    epsilon = methodParams.epsilon;
    xmin = starting_point;
    t = methodParams.startingPoint;
    nu = 0.1;
    it = 1;
    
    
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    pk = - grad;
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, pk', xmin, t, it);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        % update values
        fPrev = fCurr;
        gradOld = grad;
        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        
        % compute parameter beta
        betaFR = (grad'*grad)/(gradOld'*gradOld);
        
        % restart
        restartCoef = abs(grad'*gradOld) / (grad'*grad);
        if (restartCoef > nu)
           betaFR = 0;
        end
        
        pk = betaFR*pk - grad;
        
        it = it + 1;    
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, norm(grad));
        valuesPerIter.setStepVal(it, t);
    end

    cpuTime = toc;
    fmin = fCurr;
    valuesPerIter.trim(it);
    it = it - 1;
end
