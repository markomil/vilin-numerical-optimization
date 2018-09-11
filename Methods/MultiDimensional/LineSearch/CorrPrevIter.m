function [ outT, outX, outVal, outGr, evalNumbers ] = CorrPrevIter( functionName, params )

%%%%%%%%                Header              %%%%%%%%%%
%       This is algorithm for computing current step size parameter
%       It is computed based on the previous stepsize value 
%       as well as current function value. 
%
%       It follows simple rule:
%       if newFunVal > currFunVal
%           stepsize = coef * stepsize
%       end
%       coef < 1, usually coef = 0.5
%       
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val = vals(end); % take last (current) function value
    dir = params.dir;
    it = 1;                                         % number of iteration
    %t = params.tInitStart;                          % starting value for t
    t = params.tPrev;                               % starting value for t
    
    coef = 0.5;
    
    [val1, ~, ~] = feval(functionName, x0+t*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
          
    % process  
    while (val1 >= val)
        t = t * coef; 
        it = it + 1;
        [val1, ~, ~] = feval(functionName, x0+t*dir, [1 0 0]);
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

