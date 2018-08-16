function [ outT, outX, outVal, outGr, evalNumbers ] = Wolfe( functionName, params)
%%%%%%%%                Header              %%%%%%%%%%
%           This is Wolfe rule for 
%           inexact line search  where step size 
%           is computed by an idea
%           introduced by Nocedal and Wright 
%       
%%%%%%%%                End                 %%%%%%%%%%
	
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val0 = vals(end); % take last (current) function value
   
    gr0 = params.grad;
    dir = params.dir;
    rho = params.rho;
    sigma = params.sigma;
    ksi = params.ksi;
    tInit = params.tInitStart;
    iterNum = params.it; % number of iter of original method (outer loop)
    it = 1;              % number of iteration
    tMax = 10^(10);
    
    % This block of code determines starting value for t
    if iterNum == 1
        t = tInit;
    else
        val00 = vals(end-1); % take one before last function value
        % compute initial stepsize according to Nocedal simple rule
        t = computLineSearchStartPoint(val0, val00, gr0, dir); 
    end
       
    t1 = 0; t2 = t;                         % starting values for t0 and t1
    derPhi0 = gr0'*dir';                    % derivative of Phi(t) in  point x0
    
    % set values in point x0+t1*dir, t1 = 0
    val1 = val0;
    derPhi1 = derPhi0;
              
    while 1
        
        [val2, gr2, ~] = feval(functionName,x0+t2*dir,[1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        derPhi2 = gr2'*dir';                    % derivative of Phi(t) in current point         
               
        % check if current iterate violates sufficient decrease
        if  (val2 > val0 + derPhi0*rho*t2) || ((val2 >= val1 ) && (it > 1)) 
            % there has to be an acceptable point between t0 and t1
            % (because rho > sigma)
            [t, outVal, outGr, evalNumLocal] = nocZoom(functionName,x0,dir,val0,derPhi0,t1,t2,val1,val2,derPhi1,derPhi2,rho,sigma,ksi);
            evalNumbers = evalNumbers + evalNumLocal;
            break;
        end
        % current iterate has sufficient decrease, but are we too close?
        if (derPhi2 >= sigma*derPhi0)
            % wolfe fullfilled, quit
            t = t2;
            outVal = val2;
            outGr = gr2;
            break;
        end
        
        % update values
        t1 = t2;
        val1 = val2;
        derPhi1 = derPhi2;
        
        multCoef = 10;
        t2 = min(tMax, t2*multCoef);
        
        it = it + 1;
    end
    
    % save output values
    outX = x0 + t*dir;
    outT = t;
       
end


function [t, val1, gr1, evalNumLocal] = nocZoom(functionName,x0,dir,val0,derPhi0,tLo,tHi,valLo,valHi,derPhiLo,derPhiHi,rho,sigma,ksi)

    inter = 'cubic';
    %inter = 'quadratic';
    
    evalNumLocal = EvaluationNumbers(0, 0 ,0);

    while 1
                
        switch inter
            case 'quadratic'
                t = -0.5*(derPhiLo*(tLo^2-tHi^2) - 2*tLo*(valLo-valHi))/(valLo-valHi - derPhiLo*(tLo-tHi));
            case 'cubic'
                %qubic interpolation according to Nocedal
                if tLo < tHi 
                    t = interCubic(tLo, tHi, valLo, valHi, derPhiLo, derPhiHi);
                else
                    t = interCubic(tHi, tLo, valHi, valLo, derPhiHi, derPhiLo);
                end
            otherwise
                t = (tLo + tHi)/2;
        end;
        
        x1 = x0 + t*dir;
        [val1, gr1, ~] = feval(functionName, x1, [1 1 0]);
        evalNumLocal.incrementBy([1 1 0]);
        derPhi1 = gr1'*dir';
        
        if abs(val1 - valLo)/(1 + abs(val1)) < ksi || abs(val1 - valHi)/(1 + abs(val1)) < ksi
            return;
        end

        if ((val1 > val0 + rho*t*derPhi0) || (val1 >= valLo)) 
            % if we do not observe sufficient decrease in point t, we set
            % the maximum of the feasible interval to t
            tHi = t;
            valHi = val1;
            derPhiHi = derPhi1;
        else
            % strong wolfe fullfilled?
            if derPhi1 >= sigma*derPhi0
                return;
            end
            
            tLo = t;
            valLo = val1;
            derPhiLo = derPhi1;
        end
        
    end

end

function [t] = interCubic(t1, t2, val1, val2, der1, der2)
%   This function computes point t between points t1 and t2
%   by cubic interpolation t1 < t < t2.

    d1 = der1+der2-3*(val1-val2)/(t1-t2);
    d2 = sqrt(d1^2-der1*der2);
    tmp = t2 - (t2-t1)*(der2+d2-d1)/(der2 - der1 + 2*d2);
    
    if tmp >= 0
        d2 = sqrt(d1^2-der1*der2);
        t = t2 - (t2-t1)*(der2+d2-d1)/(der2 - der1 + 2*d2);
    else
        t = t1 - 1;
    end
    
    % if minimum is is not in the interval (t1, t2) then minimum is in t1
    if t < t1 || t > t2
        t = t1;
    end
end