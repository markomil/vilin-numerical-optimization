function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = CG_Descent( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *       CG_Descent algorithm        *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The CG_DESCENT algorithm is a nonlinear congugate-gradient
%   method for solving large-scale unconstrained minimization problem 
%   originally designed by W.W. Hager and H. Zhang. In order to achieve 
%   better numerical performances the authors suggest to use Approximate 
%   Wolfe line search which they proposed.

%   W.W. Hager H. Zhang, 
%   A new conjugate gradient method with guaranteed descent
%   and an efficient line search, 
%   SIAM J. Optim., 16(1):170–192, 2005.

%   W.W. Hager, H. Zhang,
%   Algorithm 851: "CG_Descent, a conjugate gradient method with 
%   guaranteed descent", 
%   ACM Trans. Math. Software, 32(1):113-137, 2006.

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
    tPrev = t;
    
    it = 1;
    ni = 0.01;
    
    [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, norm(grad));
    
    pk = - grad;
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1; 
    Delta = 0.7;
    Q = 0;
    C = 0;
                
    % process
    while (it < maxIter && norm(grad) > epsilon && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
    %while (it < maxIter && norm(grad, Inf) > epsilon && t*abs(grad'*pk) > 1e-20*abs(fCurr))
        
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        Q = 1 + Delta * Q;
        C = C + ((abs(fCurr) - C) / Q);
        params = LineSearchParams(methodParams, fValues, grad, pk', xmin, tPrev, it, C);
          
        [t, xmin, lineSearchEvalNumbers] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        % update values
        fPrev = fCurr;
        gradOld = grad;
        tPrev = t;
        
        [fCurr, grad, ~] = feval(functionName, xmin, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);

        % compute parameter beta
        niK = -1 / (norm(pk) * min(ni, norm(grad)));
        yk = grad - gradOld;
        py = pk' * yk;
        betaCGD = (1 / py)*(yk - 2*pk*(norm(yk)^2 / py))'*grad;
        betaCGD = max(betaCGD, niK); % Restart
        pk = betaCGD*pk - grad;
        
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
