function [ outT, outX, evalNumbers ] = Wolfe_ApproxWolfe( functionName, params )
    it = params.it;
    fvals = params.vals;
    w = params.w;
    C = params.C;
    
    if (it >= 2 && abs(fvals(end) - fvals(end - 1)) <= w*C)
        params.ksi = 1e-6 * C;
        [outT, outX, evalNumbers] = ApproxWolfe(functionName, params);
    else
        [outT, outX, evalNumbers] = Wolfe(functionName, params);
    end