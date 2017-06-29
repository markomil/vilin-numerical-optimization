function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = GoldsteinPrice( functionName, methodParams )
%%%%%%%%                Header              %%%%%%%%%%
%       This is Modified Newton method implemented by 
%       using numerical gradient and Hessian computations.
%       Parametar eta in Goldstein Price method determines 
%       wheater to use gradient or Newton search direction. 
%       Line search method is used for computing step 
%       size in every iteration.
%
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    xmin = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    eta = 0.2;                              % parametar eta initialization
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    
    [fCurr, gr, Hes] = feval(functionName, xmin, [1 1 1]);
    evalNumbers.incrementBy([1 1 1]);
    grNorm = double(norm(gr));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
                
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        % Computes dir according to the Goldstein Price rule 
        dir = (-Hes\gr)';                  % computes Newton direction
        if dir*(-gr)/(norm(dir)*norm(gr)) < eta || sum(isnan(dir)) > 0
            dir = -gr';
        end
        
        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it, gr, dir', methodParams.startingPoint);
        params = LineSearchParams(methodParams, fCurr, gr, dir, xmin, lsStartPnt);
        [t, xmin, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        it = it + 1;
            
        fPrev = fCurr;
        % compute numerical gradient and Hessian in new point
        [fCurr , gr, Hes] = feval(functionName, xmin, [1 1 1]);   
        evalNumbers.incrementBy([1 1 1]);
        grNorm = double(norm(gr));
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
    end;
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    fmin = fCurr;
    it = it - 1;
end
