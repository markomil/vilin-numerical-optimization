function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = doglegSR1( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *        Dogleg SR1 method          *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Dogleg SR1 algorithm is a trust region method for solving 
%   unconstrained minimization problem. The idea is to find the minimum 
%   of the approximation of the objective function inside some small 
%   region around the current point. Instead of using expensive Hessian
%   and it's inverse the SR1 approximation is used. Inside the trust 
%   region two directions are used steepest descent and quasi newton 
%   direction obtained by SR1 approximation.  The idea for the 
%   implementation is taken from Nocedal book 'Numerical Optimization'.

%   C.G. Broyden,
%   A class of methods for solving nonlinear simultaneous equations,
%   Mathematics of Computation, 19 (1965) 577–593.

%   J. Nocedal, S.J. Wright,
%   Numerical Optimization (2nd ed.),
%   Berlin, New York: Springer-Verlag, 2006.

%   ------------------      *******************        ------------------

    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    dim = length(x0);
    H = eye(dim);           % define initial Hes inverse approx 
    B = eye(dim);           % define initial Hes approx 
    r = 10^(-8);            % coef which determine to make update of H or not     
    
    % current xmin
    [fPrev, gr0, ~] = feval(functionName, x0, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr0));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fPrev);
    valuesPerIter.setGradientVal(it, grNorm);
    
    workPrec = methodParams.workPrec;
    fCurr = fPrev + 1;
    mPrev = fPrev;
    
    % Additional parameter initialization
    eta = 0.001; % parameter which satisfies eta in [0, 0.25)
    trustDelta = 1;
    trustDeltaMax = 10^9;
    
    
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        % computes dogleg direction with SR1 aproximation
        dir = doglegDirection(gr0, grNorm, B, H, trustDelta);
        dirNorm = norm(dir);
                
        % compute function value in current point attempt
        [fCurr, ~, ~] = feval(functionName, x0 + dir, [1 0 0]);   
        evalNumbers.incrementBy([1 0 0]);
        
        % compute function model value in current point attempt
        mCurr = quadraticModelFunction(fPrev, dir, gr0, B);
         
        % computes ratio which says whether model function is close to the original one
        rho = (fPrev - fCurr) / (mPrev - mCurr);
                        
        if rho < 0.1 
            trustDelta = 0.5 * dirNorm;
        else
            if rho > 0.75 && dirNorm > 0.8 * trustDelta
                trustDelta = min(2*trustDelta, trustDeltaMax);
            end;
        end;
        
        % update point and gradient  
        x1 = x0 + dir;
        [~, gr1, ~] = feval(functionName, x1, [0 1 0]);
        evalNumbers.incrementBy([0 1 0]);

        % update SR1 Hessian aproximation
        % compute vectors s and y
        s = (x1 - x0)';
        y = gr1 - gr0;

        % update inverse and original Hessian approximation
        if (s-H*y)'*y >= r*norm(y)*norm(s-H*y)
            H = H + (s-H*y)*(s-H*y)'/((s-H*y)'*y);
            B = B + (y-B*s)*(y-B*s)'/((y-B*s)'*s); 
        end; 
        
                
        % update current point and other data
        if rho > eta
            % update function values
            fPrev = fCurr;
            mPrev = fPrev;
            
            x0 = x1; gr0 = gr1;             % update point and gradient
            grNorm = norm(gr0);
            it = it + 1;

            valuesPerIter.setFunctionVal(it, fCurr);
            valuesPerIter.setGradientVal(it, grNorm);
                        
            fCurr = fPrev + 1;
        end;
                
    end;
    
    % determine total CPU time
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x0; 
    fmin = fPrev;
    it = it - 1;
       
end

% This function computes dogleg direction 
function outDir = doglegDirection(gr, grNorm, B, H, trustDelta)
    

    % Steepest descent direction 
    dirC = -(gr'*gr) / (gr'*B*gr) * gr';
    % Newton direction
    dirB = -(H*gr)'; 
    
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



