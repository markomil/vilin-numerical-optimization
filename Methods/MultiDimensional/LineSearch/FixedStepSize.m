function [ outT, outX, outVal, outGr, evalNumbers ] = FixedStepSize( functionName, params )

%   ------------------      *******************        ------------------
%   *                                                                   *
%   *               *************************************               *
%   *               *                                   *               *
%   *               *          Fixed step size          *               *
%   *               *                                   *               *
%   *               *************************************               *
%   *                                                                   *
%   ------------------      *******************        ------------------

%   This is fixed step size line search which means that in each iteration
%   fixed value of step size is used. 

%   ------------------      *******************        ------------------
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    dir = params.dir;
    t = params.tInitStart;                             % starting value for t
           
    % save output values
    xmin = x0 + t*dir;
    outX = xmin; outT = t;
    % compute function and gradient value in current point xmin 
    [outVal, outGr, ~] = feval(functionName, xmin, [1 1 0]);   
    evalNumbers.incrementBy([1 1 0]);
        
end

