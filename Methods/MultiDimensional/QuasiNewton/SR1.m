function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = SR1( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *             SR1 Method            *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The SR1 is quasi Newton method implemented by using numerical 
%   gradient computation and inverse Hessian aproximation. It is 
%   originally developed by Broyden. The Hessian approximation and it's 
%   inverse is obtained by symmetric rank one update (SR1) which 
%   satisfies secant equation.

%   C.G. Broyden,
%   A class of methods for solving nonlinear simultaneous equations, 
%   Mathematics of Computation, 19 (1965) 577–593.

%   ------------------      *******************        ------------------
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x1 = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    t = methodParams.startingPoint;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    dim = length(x1);
    H = eye(dim);                           % define approx matrix
    r = 10^(-8);            % coef which determine to make update of H or not                           
    
    [fCurr, gr1, ~] = feval(functionName, x1, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr1));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    % add values for plot
    if (size(x1, 2) == 2)
        valuesPerIter.setXVal(it, x1);
    end
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dir = (-H*gr1)'; % computes direction
        
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, gr1, dir, x1, t, it);
        % update values
        fPrev = fCurr; 
        x0 = x1; gr0 = gr1;
        
        % Computes x1 and step-size according to the line search method rule
        [t, x1, fCurr, gr1, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        grNorm = double(norm(gr1));
        
        % compute vectors s and y
        s = (x1 - x0)';
        y = gr1 - gr0;
        
        % update inverse Hessian approximation
        if (s-H*y)'*y >= r*norm(y)*norm(s-H*y)
            H = H + (s-H*y)*(s-H*y)'/((s-H*y)'*y); 
        end
                        
        x0 = x1; gr0 = gr1;             % update point and gradient
        it = it + 1;
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
        % add values for plot
        if (size(x1, 2) == 2)
            valuesPerIter.setXVal(it, x1);
            valuesPerIter.setDirVal(it, dir);
        end
    end
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x1; 
    fmin = fCurr;
    it = it - 1;
end
