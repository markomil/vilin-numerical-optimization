classdef Results
    % Results of optimization method
    %   fmin - min value of the function
    %   xmin - point at which minimum is found
    %   gradNorm - function gradient norm when method finished
    %   iterNum - total number of iterations
    %   cpuTime - evaluation time in seconds
    %   evaNumbers - number of evaluations of function, gradient and hessian
    %   valuesPerIteration - function, gradient, hessian value per iteration
    
    properties
        fmin, xmin, gradNorm, iterNum, cpuTime, evalNumbers, valuesPerIteration
    end
    
    methods
        function obj = Results(fmin, xmin, gradNorm, iterNum, cpuTime, evalNumbers, valuesPerIteration)
        % Constructor
            obj.fmin = fmin;
            obj.xmin = xmin;
            obj.gradNorm = gradNorm;
            obj.iterNum = iterNum;
            obj.cpuTime = cpuTime;
            obj.evalNumbers = evalNumbers;
            obj.valuesPerIteration = valuesPerIteration;
        end            
    end
    
end

