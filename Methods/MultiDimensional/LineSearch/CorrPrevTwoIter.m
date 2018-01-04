function [ outT, outX, evalNumbers ] = CorrPrevTwoIter( functionName, params )

%%%%%%%%                Header              %%%%%%%%%%
%       This is algorithm for computing current step size parameter
%       It is computed based on the previous stepsize value 
%       as well as previous and current function value. 
%       This algorithm is generalization of CorrPrevIter algorithm
%       Stepsize can be decreased as well as increased depending of the
%       information obtained from previous two iterations
%
%       It follows simple rule: 
%       if newFunVal < currFunVal && currFunVal < prevFunVal
%           stepsize = stepsize * coef1;
%       else if newFunVal > currFunVal
%               stepsize = stepsize * coef2;
%            end
%       end

%       coef1 > 1, 0 < coef2 < 1; usually coef2 = 1.5, coef2 = 0.5
%       
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val0 = vals(end); % take last (current) function value
    dir = params.dir;
    iterNum = params.it; % number of iter of original method (outer loop)
    t = params.tPrev;                              % starting value for t
    
    coef1 = 1.2;
    coef2 = 0.5;
    
    [val1, ~, ~] = feval(functionName, x0+t*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
          
    if iterNum > 1 
        % This is when we have at least two previous function values 
        val00 = vals(end-1);
        
        if val1 < val0 && val0 < val00
            t = t * coef1; % increase line search parameter
        else if val1 >= val0
                t = t * coef2; % decrease line search parameter
            end
        end
    else
        % This is if we only have one function value stored
        while (val1 >= val0)
            t = t * coef2; 
            [val1, ~, ~] = feval(functionName, x0+t*dir, [1 0 0]);
            evalNumbers.incrementBy([1 0 0]);
        end
    end
                
    % save output values
    xmin = x0 + t*dir;
    outX = xmin; outT = t;
end

