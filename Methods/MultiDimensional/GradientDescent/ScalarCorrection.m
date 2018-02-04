function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = ScalarCorrection( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *      Scalar Correction method     *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   This is gradient based two-point step size method for solving 
%   Large-scale unconstrained minimization problem originally designed 
%   by Miladinovic, Stanimirovic, Miljkovic. Similarly as for
%   Barzilai Borwein method, in order to achieve good numerical 
%   performance it is suggested to accompanied it with nonmonotone line 
%   search introduced by  L. Grippo, F. Lampariello, S. Lucidi.

%   M. Miladinovic, P. Stanimirovic, S. Miljkovic, 
%   Scalar Correction Method for Solving Large 
%   Scale Unconstrained Minimization Problems
%   J. Optim. Theory. Appl., 151 (2011) 304--320.

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
    
    [fCurr, gr0, ~] = feval(functionName, x0, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr0));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    
    gamma = 1;
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
        
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        dir = (-gamma*gr0)'; % computes direction
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, gr0, dir, x0, t, it);
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
        
        % auxilliary variable
        r = s - gamma*y;
        
        % update inverse Hessian approximation
        if y'*r > 0
            gamma = (s'*r) / (y'*r);
        else
            gamma = norm(s) / norm(y);
        end;
        
        % parameter reset in case of negative value
        if gamma < 0
            gamma = 1;
        end;
        
        x0 = x1; gr0 = gr1;             % update point and gradient
        it = it + 1;
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
    end;
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x0; 
    fmin = fCurr;
    it = it - 1;
        
end
