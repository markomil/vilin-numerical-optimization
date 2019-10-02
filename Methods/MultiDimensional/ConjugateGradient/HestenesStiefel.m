function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = HestenesStiefel( functionName, methodParams )
% Hestenes-Stiefel version of Conjugate gradient method

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *     Hestenes-Stiefel method       *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   This is nonlinear conjugate gradient method for solving 
%   large-scale unconstrained minimization problem originally designed 
%   by Hestenes and Stiefel. In order to converge Strong Wolfe line 
%   search should be applied.

%   M.R. Hestenes, E. Stiefel, 
%   Methods of Conjugate Gradients for Solving Linear Systems, 
%   J. Research Nat. Bur. Standards., 49 (1952) 409-436.

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
    nu = methodParams.nu;
    it = 1;
    
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(grad));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    % add values for plot
    if (size(xmin, 2) == 2)
        valuesPerIter.setXVal(it, xmin);
        valuesPerIter.setDirVal(it, -grad);
    end
    pk = - grad;
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    % process
    while (it < maxIter && grNorm > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, grad, pk', xmin, t, it);
        % update values
        fPrev = fCurr;
        gradOld = grad;
        
        % Computes xmin and step-size according to the line search method rule
        [t, xmin, fCurr, grad, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        grNorm = double(norm(grad));
        
        % compute parameter beta
        betaHS = (grad'*(grad-gradOld))/((grad-gradOld)'*pk);
        
        % restart
         restartCoef = abs(grad'*gradOld) / (grad'*grad);
         if (restartCoef > nu)
            betaHS = 0;
         end
        
        pk = betaHS*pk - grad;
        
        it = it + 1;
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
        % add values for plot
        if (size(xmin, 2) == 2)
            valuesPerIter.setXVal(it, xmin);
            valuesPerIter.setDirVal(it, pk);
        end
    end

    cpuTime = toc;
    fmin = fCurr;
    valuesPerIter.trim(it);
    it = it - 1;
end


