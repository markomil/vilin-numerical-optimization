function [ outT, outX, evalNumbers ] = Armijo( functionName, params )

%%%%%%%%                Header              %%%%%%%%%%
%       This is Armijo rule for inexact line search 
%       It uses quadratic and cubic interpolation
%       according to Nocedal book Numerical Optimization
%       
%%%%%%%%                End                 %%%%%%%%%%
    
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
    it = 1;                                         % number of iteration
    
    % This block of code determines starting value for t
    if iterNum == 1
        t1 = tInit;
    else
        val00 = vals(end-1); % take one before last function value
        % compute initial stepsize according to Nocedal simple rule
        t1 = computLineSearchStartPoint(val, val00, gr, dir); 
    end;
        
    [val1,~] = feval(functionName, x0+t1*dir, [1 0 0]);
    evalNumbers.incrementBy([1 0 0]);
    derPhi0 = gr'*dir';
    val2 = val1;
    t2 = t1;
          
    % process  
    while (val2 > val + rho*t2*derPhi0)
        
        if it == 1
            t2 = interQuadratic(t1, val, val1, derPhi0);
            [val2,~] = feval(functionName, x0+t2*dir, [1 0 0]);
            evalNumbers.incrementBy([1 0 0]);
        else
            [t] = interQubic(t1, t2, val, val1, val2, derPhi0);
            val1 = val2;
            t1 = t2;
            t2 = t;
            [val2,~] = feval(functionName, x0+t2*dir, [1 0 0]);
            evalNumbers.incrementBy([1 0 0]);
        end;
        it = it + 1;
         
     end; 
                
    % save output values
    xmin = x0 + t2*dir;
    outX = xmin; outT = t2;
        
end

function [t] = interQuadratic(t1, val0, val1, der0)
    t = -der0*t1^2 / (2*(val1 - val0 - der0*t1));
end

function [t] = interQubic(t1, t2, val0, val1, val2, der0)
    a = 1/((t1^2*t2^2)*(t2-t1)) * [t1^2, -t2^2] * [val2 - val0 - der0*t2; val1 - val0 - der0*t1];
    b = 1/((t1^2*t2^2)*(t2-t1)) * [-t1^3, t2^3] * [val2 - val0 - der0*t2; val1 - val0 - der0*t1];
    t = (-b + sqrt(b^2-3*a*der0)) / (3*a);
end
