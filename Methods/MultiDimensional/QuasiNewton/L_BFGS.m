function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = L_BFGS( functionName, methodParams )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *           L_BFGS Method           *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   The Limited-memory BFGS, is quasi Newton method (shortly L-BFGS) 
%   proposed by Liu and Nocedal. The method is based on the BFGS updating 
%   formula. Instead of storing fully dense matrix, only a few vectors 
%   of length n that represent the BFGS approximations are stored. 
%   Despite these optimal storage requirements, it yields good rate of 
%   convergence. In order to preserve good numerical properties Wolfe 
%   line search should be applied.

%   D.C. Liu, J. Nocedal,
%   On the limited-memory BFGS method for large scale optimization, 
%   Math. Prog., 45 (1989) 503–528.

%   ------------------      *******************        ------------------

    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x1 = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    t = methodParams.startingPoint;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
        
    [fCurr, gr1, ~] = feval(functionName, x1, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr1));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    % add values for plot
    if (size(x1, 2) == 2)
        valuesPerIter.setXVal(it, x1);
    end
    cacheSize = 5; % max number of vectors for keeping in cache
    hCoef = 1; % coef for computing initial matrix H in each iteration
    yCache = [];
    sCache = [];
    rhoCache = [];
    
    workPrec = methodParams.workPrec;
    fPrev = fCurr + 1;
        
    % process
    while (grNorm > eps && it < maxIter && abs(fPrev - fCurr)/(1 + abs(fCurr)) > workPrec)
        
        % Computes Hessian aproximation and search direction
        H = hCoef;%*eye(dim);
        dir = TwoLoopRecursion(H, gr1, sCache, yCache, rhoCache); % computes direction
                
        fValues = valuesPerIter.functionPerIteration(1:it); % take vector of function values after first 'it' iteration
        params = LineSearchParams(methodParams, fValues, gr1, dir, x1, t, it);
        % update values
        fPrev = fCurr; 
        x0 = x1; gr0 = gr1;
        
        % Computes x1 and step-size according to the line search method rule
        [t, x1, fCurr, gr1, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        grNorm = double(norm(gr1));
                              
        % compute vectors s and y and add them to cache
        s = (x1 - x0)'; sCache = addToCache(sCache, s, cacheSize);
        y = gr1 - gr0; yCache = addToCache(yCache, y, cacheSize);
        
        rho = (y'*s)^-1; rhoCache = addToCache(rhoCache, rho, cacheSize);
        hCoef = (s'*y) / (y'*y); % Nocedal, page 226
                
        it = it + 1;
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
        % add values for plot
        if (size(x1, 2) == 2)
            valuesPerIter.setXVal(it, x1);
            valuesPerIter.setDirVal(it, dir);
        end
    end
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x1; 
    fmin = fCurr;
    it = it - 1;
end
        
% auxiliary functions
function [dir] = TwoLoopRecursion(H, grad, sCache, yCache, rhoCache)
	len = size(sCache, 2); q = grad; alphas = [];

	for i = len:-1:1
	  alpha = rhoCache(:,i)*sCache(:,i)'*q;
	  alphas(i) = alpha;
	  q = q - alpha*yCache(:, i);
	end
	
	r = H*q;

	for i = 1:len
	  beta = rhoCache(:,i)*yCache(:,i)'*r;
	  r = r + sCache(:,i)*(alphas(i)-beta);
	end
	
	dir = -r';
end


function [newCache] = addToCache(cache, el, maxLength)
	len = size(cache, 2);
	if len < maxLength
	  cache(:, len+1) = el;
	else
	  cache(:, 1) = [];
	  cache(:, maxLength) = el;
	end
	newCache = cache;
end