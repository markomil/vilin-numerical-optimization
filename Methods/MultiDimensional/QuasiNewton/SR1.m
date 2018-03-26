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
    x0 = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    t = methodParams.startingPoint;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    dim = length(x0);
    H = eye(dim);                           % define approx matrix
    r = 10^(-8);            % coef which determine to make update of H or not                           
    
    [fCurr, gr0, ~] = feval(functionName, x0, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr0));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
    
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        % Computes xmin according to the method rule 
        dir = (-H*gr0)';                    % computes direction
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, gr0, dir, x0, t, it);
        [t, x1, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
            
        % update function value
        fPrev = fCurr;
        % compute numerical gradient in new point
        [fCurr, gr1, ~] = feval(functionName, x1, [1 1 0]);   
        evalNumbers.incrementBy([1 1 0]);
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
    end
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x0; 
    fmin = fCurr;
    it = it - 1;
end
