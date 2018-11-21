function [ outT, outX, outVal, outGr, evalNumbers ] = Backtracking( functionName, params )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *     Backtracking line search      *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   Backtracking line search is a line search procedure for computing 
%   step-size parameter t such that it satisfies so called Armijo rule.
%   This method ensures a sufficient decrease of the goal function.   
%   For finding such a parameter standard backtracing technique is used. 
%   This idea is originally proposed by L. Armijo

%   L. Armijo,
%   Minimization of functions having Lipschitz first partial derivatives,
%   Pac. J. Math, 6 (1966) 1-3.

%   ------------------      *******************        ------------------
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val = vals(end); % take last (current) function value
    gr = params.grad;
    dir = params.dir;
    rho = params.rho;
    beta = params.beta;
    tInit = params.tInitStart;
    iterNum = params.it; % number of iter of original method (outer loop)
    it = 1;                                         % number of iteration
    
    % This block of code determines starting value for t
    if iterNum == 1
        t = tInit;
    else
        val00 = vals(end-1); % take one before last function value
        % compute initial stepsize according to Nocedal simple rule
        t = computLineSearchStartPoint(val, val00, gr, dir); 
    end;
        
    [val1,~] = feval(functionName, x0+t*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
    derPhi0 = gr'*dir';
      
    % process  
    while (val1 > val + rho*t*derPhi0)
        
        t = t * beta; 
        it = it + 1;
        [val1,~] = feval(functionName, x0+t*dir, [1 0 0]);
        evalNumbers.incrementBy([1 0 0]);
    end; 
                
    % save output values
    xmin = x0 + t*dir;
    outX = xmin; outT = t;
    outVal = val1;
    % compute gradient in current point xmin 
    [~, outGr, ~] = feval(functionName, xmin, [0 1 0]);   
    evalNumbers.incrementBy([0 1 0]);
    
end