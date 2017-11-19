function [ fmin, xmin, it, cpuTime, evalNumbers, valuesPerIter ] = L_BFGS( functionName, methodParams )
%%%%%%%%                Header              %%%%%%%%%%
%       This is "Limited memory" version of BFGS method.
%	This is "TwoLoopRecursion" version of the method.
%
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = methodParams.starting_point;
    maxIter = methodParams.max_iteration_no;
    valuesPerIter = PerIteration(maxIter);
    eps = methodParams.epsilon;
    tic;                                    % to compute CPU time
    it = 1;                                 % number of iteration
    dim = length(x0);
    
    [fCurr, gr0, ~] = feval(functionName, x0, [1 1 0]);
    evalNumbers.incrementBy([1 1 0]);
    grNorm = double(norm(gr0));
    % Added values for first iteration in graphic
    valuesPerIter.setFunctionVal(it, fCurr);
    valuesPerIter.setGradientVal(it, grNorm);
    
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
        H = hCoef; %*eye(dim);
        dir = TwoLoopRecursion(H, gr0, sCache, yCache, rhoCache); % computes direction
        % Computes xmin according to the method rule
        lsStartPnt = computLineSearchStartPoint(fCurr, fPrev, it, gr0, dir', methodParams.startingPoint);
        params = LineSearchParams(methodParams, fCurr, gr0, dir, x0, lsStartPnt);
        [t, x1, lineSearchEvalNumbers ] = feval(methodParams.lineSearchMethod, functionName, params);
        evalNumbers = evalNumbers + lineSearchEvalNumbers;
        
        % update function value
        fPrev = fCurr;
        % compute numerical gradient in new point
        [fCurr, gr1] = feval(functionName, x1, [1 1 0]);   
        evalNumbers.incrementBy([1 1 0]);
        grNorm = double(norm(gr1));
        
        % compute vectors s and y and add them to cache
        s = (x1 - x0)'; sCache = addToCache(sCache, s, cacheSize);
        y = gr1 - gr0; yCache = addToCache(yCache, y, cacheSize);
        
        rho = (y'*s)^-1; rhoCache = addToCache(rhoCache, rho, cacheSize);
        hCoef = (s'*y) / (y'*y); % Nocedal, page 226
        
        x0 = x1; gr0 = gr1;             % update point and gradient
        it = it + 1;
        
        valuesPerIter.setFunctionVal(it, fCurr);
        valuesPerIter.setGradientVal(it, grNorm);
        valuesPerIter.setStepVal(it, t);
    end;
    
    cpuTime = toc;
    valuesPerIter.trim(it);
    xmin = x0; 
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