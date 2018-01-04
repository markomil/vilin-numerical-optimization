function [ outT, outX, evalNumbers ] = Goldstein( functionName, params )
%%%%%%%%          Header              %%%%%%%%%%
%       This is Goldstein rule for 
%           inexact line search 
%       
%%%%%%%%           End                %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val = vals(end); % take last (current) function value
    gr = params.grad;
    dir = params.dir;
    rho = params.rho;
    tInit = params.tInitStart;
    iterNum = params.it; % number of iter of original method (outer loop)
    it = 1;                                 % number of iteration

    % This block of code determines starting value for t
    if iterNum == 1
        t = tInit;
    else
        val00 = vals(end-1); % take one before last function value
        % compute initial stepsize according to Nocedal simple rule
        t = computLineSearchStartPoint(val, val00, gr, dir); 
    end;
    
    gamma = 1.1;                            % set value for gamma
    t1 = 0; t2 = Inf;                       % starting values for t1 and t2
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

