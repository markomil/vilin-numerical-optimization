function [ outT, outX, evalNumbers ] = ApproxWolfe1( functionName, params )    
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
    [aj, bj, evalNumbersB, valAj, derAj, valBj, derBj] = bracket(c, val0, derPhi0, functionName, x0, dir, 5, theta, eps);
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
                    
        [a, b, evalNumbersS2, valA, derA, valB, derB] = secant2(aj, bj, valAj, derAj, valBj, derBj, functionName, val0, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbersS2; 
            
        if b-a > gamma * (bj - aj)
            c = (a + b) / 2;
            [a, b, evalNumbersU, valA, derA, valB, derB] = update(a, b, c, valA, derA, valB, derB, val0, functionName, x0, dir, theta, eps);
            evalNumbers = evalNumbers + evalNumbersU;
        end
            
        aj = a;
        bj = b;
        valAj = valA;
        derAj = derA;
        valBj = valB;
        derBj = derB;
        
        
        c = min(tMax, c);
        
        it = it + 1;
    end
    
    % save output values
    outX = x0 + t*dir;
    outT = t;
       
end

function [a_, b_, evalNumbers, valA_, derA_, valB_, derB_] = update3(a, b, valA, derPhiA, valB, derPhiB, phi0, functionName, x0, dir, theta, eps)

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
            valA_ = valA;
            derA_ = derPhiA;
            valB_ = phiD;
            derB_ = derPhiD;
            break;
        end

        % U3b
        if derPhiD <  0 && phiD <= phi0 + eps
            a_ = d;
            valA_ = phiD;
            derA_ = derPhiD;
            valB_ = valB;
            derB_ = derPhiB;
        end
        
        % U3c
        if derPhiD < 0 && phiD > phi0 + eps
            b_ = d;
            valA_ = valA;
            derA_ = derPhiA;
            valB_ = phiD;
            derB_ = derPhiD;
        end
    end
end

function [a_, b_, evalNumbers, valA_, derA_, valB_, derB_] = update(a, b, c, valA, derA, valB, derB, phi0, functionName, x0, dir, theta, eps)
    evalNumbers = EvaluationNumbers(0,0,0);

    % U0
    if c <= a || c >= b
        a_ = a;
        b_ = b;
        valA_ = valA;
        derA_ = derA;
        valB_ = valB;
        derB_ = derB;
        return;
    end
    
    [phiC, derPhiC, ~] = feval(functionName, x0+c*dir, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    derPhiC =  derPhiC'*dir';
    
    % U1
    if derPhiC >= 0
        a_ = a;
        b_ = c;
        valA_ = valA;
        derA_ = derA;
        valB_ = phiC;
        derB_ = derPhiC;
        return;
    end
    
    % U2
    if derPhiC < 0 && phiC <= phi0 + eps
        a_ = c;
        b_ = b;
        valA_ = phiC;
        derA_ = derPhiC;
        valB_ = valB;
        derB_ = derB;
        return;
    end
    
    % U3
    if derPhiC < 0 && phiC > phi0 + eps
        [a_, b_, evalNumbers3, valA_, derA_, valB_, derB_] = update3(a, c, valA, derA, phiC, derPhiC, phi0, functionName, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbers3;
        return;
    end
end

function [a0, b0, evalNumbers, valA_, derA_, valB_, derB_] = bracket(c, phi0, derPhi0, functionName, x0, dir, range_expansion, theta, eps)
    cj = c;
    ci = 0;
    evalNumbers = EvaluationNumbers(0,0,0);
    
    valCI = phi0; derCI = derPhi0;
        
    while 1
        [phiJ, derPhiJ, ~] = feval(functionName, x0+cj*dir, [1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        derPhiJ =  derPhiJ'*dir';
                       
        if phiJ <= phi0 + eps
            ci = cj;
            valCI = phiJ;
            derCI = derPhiJ;
        end
        
        if derPhiJ >= 0
            b0 = cj;
            a0 = ci;
            valA_ = valCI;
            derA_ = derCI;
            valB_ = phiJ;
            derB_ = derPhiJ;
            break;
        end
        
        if derPhiJ < 0 && phiJ > phi0 + eps
            [a0, b0, evalNumbers3, valA_, derA_, valB_, derB_] = update3(0, cj, phi0, derPhi0, phiJ, derPhiJ, phi0, functionName, x0, dir, theta, eps);
            evalNumbers = evalNumbers + evalNumbers3;
            break;
        end
        
        cj = range_expansion * cj;
    end
end

function [c] = secant(a, b, derPhiA, derPhiB)
       
    d = (derPhiB - derPhiA);
    if d == 0 || isnan(d) || d == Inf || d == -Inf
        d = 1e-16;
    end
    n = (a*derPhiB - b*derPhiA);
    c = n / d;
end

function [a_, b_, evalNumbers, valA_, derA_, valB_, derB_] = secant2(a, b, valA, derA, valB, derB, functionName, phi0, x0, dir, theta, eps)
    c = secant(a, b, derA, derB);
    evalNumbers = EvaluationNumbers(0,0,0);
    
    [A, B, evalNumbersU, valA_, derA_, valB_, derB_] = update(a, b, c, valA, derA, valB, derB, phi0, functionName, x0, dir, theta, eps);
    evalNumbers = evalNumbers + evalNumbersU;
    
    if c == B 
        s = b;
        S = B;
        der_s = derB;
        der_S = derB_;
    end
    
    if c == A
        s = a;
        S = A;
        der_s = derA;
        der_S = derA_;
    end
    
    if c == A || c == B        
        c_ = secant(s, S, der_s , der_S);
        [a_, b_, evalNumbersU, valA_, derA_, valB_, derB_] = update(A, B, c_, valA_, derA_, valB_, derB_, phi0, functionName, x0, dir, theta, eps);
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