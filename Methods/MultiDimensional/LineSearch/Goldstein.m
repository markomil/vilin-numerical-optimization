function [ outT, outX, evalNumbers ] = Goldstein( functionName, params )
%%%%%%%%          Header              %%%%%%%%%%
%       This is Goldstein rule for 
%           inexact line search 
%       
%%%%%%%%           End                %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    val = params.val;
    gr = params.grad;
    dir = params.dir;
    rho = params.rho;
    it = 1;                                 % number of iteration

    gamma = 1.1;                            % set value for gamma
    t1 = 0; t2 = Inf;                       % starting values for t1 and t2
    t = params.tStart;                             % starting value for t
    [val1,~] = feval(functionName, x0 + t*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
      
    % process  
    while (val1 > val + rho*t*gr'*dir' || val1 < val + (1-rho)*t*gr'*dir')
        
        if (val1 > val + rho*t*gr'*dir')
            t2 = t; t = (t1 + t2)/2;
        else
            if (val1 < val + (1-rho)*t*gr'*dir')
                t1 = t;
                if (t2 < Inf)
                    t = (t1 + t2)/2;
                else
                    t = t * gamma;
                end
            end
        end
        it = it + 1;
        [val1, ~] = feval(functionName, x0+t*dir, [1 0 0]);
        evalNumbers.incrementBy([1 0 0]);
    end; 
                
    % save and print output values
    xmin = x0 + t*dir;
    outX = xmin; outT = t;
        
end

