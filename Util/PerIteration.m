classdef PerIteration < handle
    %Holds info about function, gradient, hessian value per iteration
    
    properties
        iterations, functionPerIteration, gradientPerIteration, hessianPerIteration, stepPerIteration, xPerIteration, directionPerIteration
    end
    
    methods
        function obj = PerIteration(MaxIterationNumber)
            %Constructor
            obj.iterations = linspace(0, MaxIterationNumber, MaxIterationNumber + 1);
            obj.functionPerIteration = zeros(1, MaxIterationNumber + 1);
            obj.gradientPerIteration = zeros(1, MaxIterationNumber + 1);
            obj.hessianPerIteration = zeros(1, MaxIterationNumber + 1);
            obj.stepPerIteration = zeros(1, MaxIterationNumber + 1);
            obj.xPerIteration = zeros(2, MaxIterationNumber + 1);
            obj.directionPerIteration = zeros(2, MaxIterationNumber + 1);
        end
        
        function setFunctionVal(obj, iter, val)
            %Sets function value in iter iteration
            obj.functionPerIteration(iter) = val;
        end
        
        function setGradientVal(obj, iter, val)
            %Sets gradient value in iter iteration
            obj.gradientPerIteration(iter) = val;
        end
        
        function setHessianVal(obj, iter, val)
            %Sets hessian value in iter iteration
            obj.hessianPerIteration(iter) = val;
        end
        
        function setStepVal(obj, iter, val)
            %Sets hessian value in iter iteration
            obj.stepPerIteration(iter) = val;
        end
        
        function setXVal(obj, iter, val)
            %Sets x value in iter iteration
            obj.xPerIteration(:, iter) = val;
        end
        
        function setDirVal(obj, iter, val)
            %Sets direction value in iter iteration
            obj.directionPerIteration(:, iter) = val;
        end
        
        function trim(obj, iterNum)
            %Keeps only values for first iterNum itartions
            obj.iterations = obj.iterations(1:iterNum);
            obj.functionPerIteration = obj.functionPerIteration(1:iterNum);
            obj.gradientPerIteration = obj.gradientPerIteration(1:iterNum);
            obj.hessianPerIteration = obj.hessianPerIteration(1:iterNum);
            obj.stepPerIteration = obj.stepPerIteration(1:iterNum);
            obj.xPerIteration = obj.xPerIteration(:, 1:iterNum);
            obj.directionPerIteration = obj.directionPerIteration(:, 1:iterNum);
        end
        
    end
    
end

