function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = DFP( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *             DFP Method            *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The DFP is quasi Newton method proposed originally by Davidon and 
%   later developed by Fletcher and Powell. This is one of the most
%   popular rank two update quasi Newton formula. In order to maintain
%   good search directions the Wolfe or strong Wolfe line search should 
%   be applied.

%   W.C. Davidon,
%   Variable metric method for minimization, 
%   SIAM J. Optim. 1 (1991) 1–17.

%   R. Fletcher, M.J.D. Powell
%   A rapid convergent descent method for minimization, 
%   Computer Journal 6 (1963) 163–168.

%   ------------------      *******************        ------------------

    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    t = methodParams.startingPoint;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    dim = length(x0);
    H = eye(dim);                           % define approx matrix
    
    [fCurr, gr0, ~] = feval(functionName,x0,[1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr0));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
                             
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dir = (-H*gr0)';                    % computes direction
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, gr0, dir, x0, t, it);
        % Computes xmin according to the method rule 
        [t, x1, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
            
        % update function value
        fPrev = fCurr;
        % compute numerical gradient in new point
        [fCurr, gr1] = feval(functionName, x1, [1 1 0]);   
        evalNumbers.incrementBy([1 1 0]);
        grNorm = double(norm(gr1));
        
        % compute vectors s and y
        s = (x1 - x0)';
        y = gr1 - gr0;
        
        % update inverse Hessian approximation
        H = H + s*s'/(s'*s) - (H*y)*(y'*H)/(y'*H*y);      
        
        x0 = x1; gr0 = gr1;             % update point and gradient
        it = it + 1;
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
    end
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x0; 
    fmin = fCurr;
    it = it - 1;
end
