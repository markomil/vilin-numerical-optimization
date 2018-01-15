function [ outT, outX, evalNumbers ] = ApproxWolfe( functionName, params )    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    vals = params.vals;
    val0 = vals(end); % take last (current) function value
    
    gr0 = params.grad;
    dir = params.dir;
    rho = params.rho; % delta in paper
    theta = params.theta;
    gamma = params.gamma;
    sigma = params.sigma;
    tInit = params.tInitStart;
    iterNum = params.it; % number of iter of original method (outer loop)
    it = 1;                               % number of iteration
    tMax = 10^(10);
    eps = params.ksi;
       
    derPhi0 = gr0'*dir';                    % derivative of Phi(t) in  point x0
    
    c = initial(x0, gr0, val0, iterNum, tInit);
    [aj, bj, evalNumbersB] = bracket(c, val0, functionName, x0, dir, 5, theta, eps);
    evalNumbers = evalNumbers + evalNumbersB;
              
    while 1
        [val2, gr2, ~] = feval(functionName,x0+c*dir,[1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        derPhi2 = gr2'*dir';                    % derivative of Phi(t) in current point         
        
        if (rho*derPhi0*c >= (val2 - val0) && derPhi2 >= sigma*derPhi0) || ... 
           (((2*rho - 1)*derPhi0 >= derPhi2 && derPhi2 >= sigma*derPhi0) || val2 <= val0 + eps)
            t = c;
            break;
        end
                    
        [a, b, evalNumbersS2] = secant2(aj, bj, functionName, val0, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbersS2; 
            
        if b-a > gamma * (bj - aj)
            c = (a + b) / 2;
            [a, b, evalNumbersU] = update(a, b, c, val0, functionName, x0, dir, theta, eps);
            evalNumbers = evalNumbers + evalNumbersU;
        end
            
        aj = a;
        bj = b;    
        
        c = min(tMax, c);
        
        it = it + 1;
    end
    
    % save output values
    outX = x0 + t*dir;
    outT = t;
       
end

function [a_, b_, evalNumbers] = update3(a, b, phi0, functionName, x0, dir, theta, eps)
    evalNumbers = EvaluationNumbers(0,0,0);
    a_ = a;
    b_ = b;
    
    while 1
        d = (1 - theta) * a_ + theta * b_;
        [phiD, derPhiD, ~] = feval(functionName, x0 + d*dir, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        derPhiD = derPhiD' * dir';

        % U3a
        if derPhiD >= 0
            b_ = d;
            break;
        end

        % U3b
        if derPhiD <  0 && phiD <= phi0 + eps
            a_ = d;
        end
        
        % U3c
        if derPhiD < 0 && phiD > phi0 + eps
            b_ = d;
        end
    end
end

function [a_, b_, evalNumbers] = update(a, b, c, phi0, functionName, x0, dir, theta, eps)
    evalNumbers = EvaluationNumbers(0,0,0);

    % U0
    if c <= a || c >= b
        a_ = a;
        b_ = b;
        return;
    end
    
    [phiC, derPhiC, ~] = feval(functionName, x0+c*dir, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    derPhiC =  derPhiC'*dir';
    
    % U1
    if derPhiC >= 0
        a_ = a;
        b_ = c;
        return;
    end
    
    % U2
    if derPhiC < 0 && phiC <= phi0 + eps
        a_ = c;
        b_ = b;
        return;
    end
    
    % U3
    if derPhiC < 0 && phiC > phi0 + eps
        [a_, b_, evalNumbers3] = update3(a, c, phi0, functionName, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbers3;
        return;
    end
end

function [a0, b0, evalNumbers] = bracket(c, phi0, functionName, x0, dir, range_expansion, theta, eps)
    cj = c;
    ci = 0;
    evalNumbers = EvaluationNumbers(0,0,0);
    
    while 1
        [phiJ, derPhiJ, ~] = feval(functionName, x0+cj*dir, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        derPhiJ =  derPhiJ'*dir';
        
        if phiJ <= phi0 + eps
            ci = cj;
        end
        
        if derPhiJ >= 0
            b0 = cj;
            a0 = ci;
            break;
        end
        
        if derPhiJ < 0 && phiJ > phi0 + eps
            [a0, b0, evalNumbers3] = update3(0, cj, phi0, functionName, x0, dir, theta, eps);
            evalNumbers = evalNumbers + evalNumbers3;
            break;
        end
        
        cj = range_expansion * cj;
    end
end

function [c, evalNumbers] = secant(a, b, functionName, x0, dir)
    evalNumbers = EvaluationNumbers(0,0,0);
    
    [~, derPhiA, ~] = feval(functionName, x0+a*dir, [0 1 0]);
    evalNumbers.incrementBy([0 1 0]);
    derPhiA =  derPhiA'*dir';
    
    [~, derPhiB, ~] = feval(functionName, x0+b*dir, [0 1 0]);
    evalNumbers.incrementBy([0 1 0]);
    derPhiB =  derPhiB'*dir';
    
    d = (derPhiB - derPhiA);
    if d == 0 || isnan(d) || d == Inf || d == -Inf
        d = 1e-16;
    end
    n = (a*derPhiB - b*derPhiA);
    c = n / d;
end

function [a_, b_, evalNumbers] = secant2(a, b, functionName, phi0, x0, dir, theta, eps)
    c = secant(a, b, functionName, x0, dir);
    evalNumbers = EvaluationNumbers(0,0,0);
    
    [A, B, evalNumbersU] = update(a, b, c, phi0, functionName, x0, dir, theta, eps);
    evalNumbers = evalNumbers + evalNumbersU;
    
    if c == B 
        s = b;
        S = B;
    end
    
    if c == A
        s = a;
        S = A;
    end
    
    if c == A || c == B
        [c_, evalNumbersS] = secant(s, S, functionName, x0, dir);
        evalNumbers = evalNumbers + evalNumbersS;
        [a_, b_, evalNumbersU] = update(A, B, c_, phi0, functionName, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbersU;
    else
        a_ = A;
        b_ = B;
    end
end

function c = initial(x0, gr0, val0, k, cOld)

ksi0 = 0.01;
ksi2 = 2;

if k == 1 % we count iterations from 1
    if x0 ~= zeros(1, length(x0))
        c = ksi0 * norm(x0, Inf) / norm(gr0, Inf);
        return;
    end

    if val0 ~= 0
        ngr0 = norm(gr0);
        c = ksi0 * abs(val0) / ngr0^2;
        return
    end

    c = 1;
    return;
end

c = ksi2 * cOld;
return;
end