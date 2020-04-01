function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = BarzilaiBorwein( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *      Barzilai Borwein method      *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Barzilai Borwein method is gradient based two point step-size 
%   method, originally developed by J. Barzilai and J.M. Borwein.
%   The idea for computing search direction comes from the 
%   secant equation. In order to achieve good numerical performance
%   Raydan suggested to accompanied it with nonmonotone line search 
%   introduced by  L. Grippo, F. Lampariello, S. Lucidi.

%   J. Barzilai and J.M. Borwein,
%   Two point step size gradient method, 
%   IMA J. Numer. Anal., 8 (1988) 141–148.

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
    gamma = 1;
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
        
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dir = (-gamma*gr1)'; % computes search direction
                
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, gr1, dir, x1, t, it);
        % update values
        fPrev = fCurr;
        gr0 = gr1;
        x0 = x1;
        
        % Computes x1 and step-size according to the line search method rule
        [t, x1, fCurr, gr1, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        grNorm = double(norm(gr1));
               
        % compute vectors s and y
        s = (x1 - x0)';
        y = gr1 - gr0;
        
        % update inverse Hessian approximation
        %gamma = (s'*y) / (s'*s); dual formula
        gamma = (s'*y) / (y'*y);
        if gamma < 0
            gamma = 1;
        end
        
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
