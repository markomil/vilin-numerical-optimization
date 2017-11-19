function [ outT, outX, evalNumbers ] = ApproxWolfe( functionName, params )    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    val0 = params.val;
    gr0 = params.grad;
    dir = params.dir;
    rho = params.rho; % delta in paper
    theta = params.theta;
    eps = params.eps;
    gama = params.gama;
    sigma = params.sigma;
    it = 1;                               % number of iteration
    tMax = 10^(10);
    t = params.tStart;                      % starting value for t
    derPhi0 = gr0'*dir';                    % derivative of Phi(t) in  point x0
    
    c = initial(x0, gr0, val0); %t);
    [aj, bj, evalNumbersB] = bracket(c, val0, functionName, x0, dir, theta, 5, eps);
    evalNumbers = evalNumbers + evalNumbersB;
              
    while 1
        [val2, gr2, ~] = feval(functionName,x0+c*dir,[1 1 0]);
        evalNumbers.incrementBy([1 1 0]);
        derPhi2 = gr2'*dir';                    % derivative of Phi(t) in current point         
               
        % check if current iterate violates sufficient decrease
        if  (((val2 > val0 + derPhi0*rho*c) && (derPhi2 >= sigma*derPhi0)) ... 
                || ((2*rho - 1)*derPhi0 < derPhi2 && val2 <= val0 + eps)) && it > 1
            t = c;
            break;
        end
                    
        [a, b, evalNumbersS2] = secant2(aj, bj, functionName, val0, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbersS2; 
            
        if b-a > gama * (bj - aj)
            c = (a + b) / 2;
            [a, b, evalNumbersU] = update(a, b, c, val0, functionName, x0, dir, theta, eps);
            evalNumbers = evalNumbers + evalNumbersU;
        end
            
        aj = a;
        bj = b;    
        
        multCoef = 10;
        c = min(tMax, c*multCoef);
        
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
        [phiD, derPhiD, ~] = eval(functionName, x0 + d*dir, [1 1 0]);
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
    if c < a || c > b
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

function [a0, b0, evalNumbers] = bracket(c, phi0, functionName, x0, dir, rho, theta, eps)
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
        
        cj = rho * cj;
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
    if d == 0
        d = 1e-8;
    end
    
    c = (a*derPhiB - b*derPhiA) / d;
end

function [a_, b_, evalNumbers] = secant2(a, b, functionName, phi0, x0, dir, theta, eps)
    c = secant(a, b, functionName, x0, dir);
    evalNumbers = EvaluationNumbers(0,0,0);
    
    [A, B, evalNumbersU] = update(a, b, c, phi0, functionName, x0, dir, theta, eps);
    evalNumbers = evalNumbers + evalNumbersU;
    
    c_ = 0;
    if c == B
        [c_, evalNumbersS] = secant(b, B, functionName, x0, dir);
        evalNumbers = evalNumbers + evalNumbersS;
    end
    
    if c == A
        [c_, evalNumbersS] = secant(a, A, functionName, x0, dir);
        evalNumbers = evalNumbers + evalNumbersS;
    end
    
    if c == A || c == B
        [a_, b_, evalNumbersU] = update(A, B, c_, phi0, functionName, x0, dir, theta, eps);
        evalNumbers = evalNumbers + evalNumbersU;
    else
        a_ = A;
        b_ = B;
    end
end

function c = initial(x0, gr0, val0)
c=1;
return
ksi0 = 0.01;

if x0 ~= zeros(1, length(x0))
    c = ksi0 * norm(x0, Inf) / norm(gr0, Inf);
    return;
end

if val0 ~= 0
    c = ksi0 * abs(val0) / norm(gr0);
    return
end

c = 1;
return;
end