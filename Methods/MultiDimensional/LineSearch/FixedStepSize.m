function [ outT, outX, evalNumbers ] = FixedStepSize( ~, params )

%%%%%%%%                Header              %%%%%%%%%%
%       This is Line Search with fixed step-size 
%       
%%%%%%%%                End                 %%%%%%%%%%
    
    % set initial values
    evalNumbers = EvaluationNumbers(0,0,0);
    x0 = params.startingPoint;
    dir = params.dir;
    t = params.tInitStart;                             % starting value for t
           
    % save output values
    xmin = x0 + t*dir;
    outX = xmin; outT = t;
        
end

