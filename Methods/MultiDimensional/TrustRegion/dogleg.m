function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = dogleg( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *           Dogleg method           *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Dogleg algorithm is a trust region method for solving 
%   unconstrained minimization problem. It is originally invented by 
%   Powell. The idea is to find the minimum of the approximation 
%   of the objective function inside some small region around the 
%   current point. Two directions are used steepest descent as well as
%   Newton direction.

%   M.J.D. Powell,
%   A new algorithm for unconstrained optimization 
%   in: J.B. Rosen, O.L. Mangasarian and K. Ritter, eds., 
%   Nonlinear Programming (1970) 31–66.

%   ------------------      *******************        ------------------

    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    xmin = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    %xmin = x0;                              % current xmin
    [fPrev, gr, Hes] = feval(functionName, xmin, [1 1 1]);
    evalNumbers.incrementBy([1 1 1]);
    grNorm = double(norm(gr));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fPrev);
    valuesPerIter.setGradientVal(it, grNorm);
    
    workPrec = methodParams.workPrec;
    fCurr = fPrev + 1;
    mPrev = fPrev;
    
    % Additional parameter initialization
    eta = 0.001; % parameter which satisfies eta in [0, 0.25)
    trustDelta = 5;
    trustDeltaMax = 10^9;
    
    
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)

        % computes dogleg direction 
        dir = doglegDirection(gr, grNorm, Hes, trustDelta);
        dirNorm = norm(dir);
                
        % compute function value in current point attempt
        [fCurr, ~, ~] = feval(functionName, xmin + dir, [1 0 0]);   
        evalNumbers.incrementBy([1 0 0]);
        
        % compute function model value in current point attempt
        mCurr = quadraticModelFunction(fPrev, dir, gr, Hes);
         
        % computes ratio which says whether model function is close to the original one
        rho = (fPrev - fCurr) / (mPrev - mCurr);
                        
        if rho < 0.1 %&& abs(fPrev - fCurr) > delta
            trustDelta = 0.25 * dirNorm;
        else
            if rho > 0.75 && abs(dirNorm - trustDelta) < eps
                trustDelta = min(2*trustDelta, trustDeltaMax);
            end
        end
        
        % update current point and all data
        if rho > eta
                        
            % update function values
            fPrev = fCurr;
            mPrev = fPrev;
            
            % update point gradient and Hessian 
            xmin = xmin + dir;
            [~, gr, Hes] = feval(functionName, xmin, [0 1 1]);
            evalNumbers.incrementBy([0 1 1]);
            grNorm = norm(gr);
                                                
            it = it + 1;

            valuesPerIter.setFunctionVal(it, fCurr);
            valuesPerIter.setGradientVal(it, grNorm);
                        
            fCurr = fPrev + 1;
        end
                
    end
    
    % determine total CPU time
    cpuTime = toc;
    
    valuesPerIter.trim(it);
    fmin = fPrev;
    it = it - 1;
        
end

% This function computes dogleg direction 
function outDir = doglegDirection(gr, grNorm, Hes, trustDelta)
    

    % Steepest descent direction 
    dirC = -(gr'*gr) / (gr'*Hes*gr) * gr';
    % Newton direction
    dirB = -(Hes\gr)'; 
    
    if norm(dirB) <= trustDelta
        outDir = dirB;
    else
        
        if norm(dirC) > trustDelta
            outDir = -trustDelta*gr'/grNorm;
        else
            tau = computeTau(dirB, dirC, trustDelta);
            outDir = dirC + (tau-1)*(dirB-dirC);
        end
    end
  
end

function outTau = computeTau(dirB, dirC, trustDelta)
    coef = zeros(1, 3);
    coef(1) = (dirB-dirC)*(dirB-dirC)';
    coef(2) = 2*(2*dirC-dirB)*(dirB-dirC)';
    coef(3) = (2*dirC-dirB)*(2*dirC-dirB)' - trustDelta^2;
    
    r = roots(coef);
    outTau = r(logical(r > 1));
end

% This is quadratic model function for objective function given by functionName
function outVal = quadraticModelFunction(fCurr, dir, gr, Hes)
    
    outVal = fCurr + dir*gr + 0.5*dir*Hes*dir';
    
end



