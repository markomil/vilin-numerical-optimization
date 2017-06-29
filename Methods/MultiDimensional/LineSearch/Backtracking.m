function [ outT, outX, evalNumbers ] = Backtracking( functionName, params )

%%%%%%%%                Header              %%%%%%%%%%
%       This is Backtracking algorithm for 
%       satisfying Armijo rule for inexact line search 
%       
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    val = params.val;
    gr = params.grad;
    dir = params.dir;
    rho = params.rho;
    beta = params.beta;
    it = 1;                                         % number of iteration
    t = params.tStart;                              % starting value for t
    
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
end