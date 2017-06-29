classdef EvaluationNumbers < handle
%Represents number of evaluations of function, gradient and hessian    
    properties
        functionEvalNo, gradientEvalNo, hessianEvalNo
    end
    
    methods	
        function obj = EvaluationNumbers(FunctionEvalNo, GradientEvalNo, HessianEvalNo)
        %Constructor
            obj.functionEvalNo = FunctionEvalNo;
            obj.gradientEvalNo = GradientEvalNo;
            obj.hessianEvalNo = HessianEvalNo;
        end
        
        function sum = plus(EN1, EN2)
        %Sum of two EvaluationNumbers objects
            sum = EvaluationNumbers(EN1.functionEvalNo + EN2.functionEvalNo, EN1.gradientEvalNo + EN2.gradientEvalNo, EN1.hessianEvalNo + EN2.hessianEvalNo);
        end
        
        function diff = minus(EN1, EN2)
        %Substraction of two EvaluationNumbers objects
            diff = EvaluationNumbers(EN1.functionEvalNo - EN2.functionEvalNo, EN1.gradientEvalNo - EN2.gradientEvalNo, EN1.hessianEvalNo - EN2.hessianEvalNo);
        end
        
        function e = eq(EN1, EN2)
        %Equality check of two EvaluationNumbers objects
            e = EvaluationNumbers(EN1.functionEvalNo == EN2.functionEvalNo, EN1.gradientEvalNo == EN2.gradientEvalNo, EN1.hessianEvalNo == EN2.hessianEvalNo);
        end
        
        function incrementFunctionEvalNo(obj)
        %Increments functionEvalNo by one
            aux = obj.functionEvalNo;
            obj.functionEvalNo  = aux + 1;
        end
        
        function decrementFunctionEvalNo(obj)
        %Decrements functionEvalNo by one
            aux = obj.functionEvalNo;
            obj.functionEvalNo  = aux - 1;
        end
        
        function incrementGradientEvalNo(obj)
        %Increment gradientEvalNo by one
            aux = obj.gradientEvalNo;
            obj.gradientEvalNo  = aux + 1;
        end
        
        function decrementGradientEvalNo(obj)
        %Decrements gradientEvalNo by one
            aux = obj.gradientEvalNo;
            obj.gradientEvalNo  = aux - 1;
        end
        
        function incrementHessianEvalNo(obj)
        %Increments hessianEvalNo by one
            aux = obj.hessianEvalNo;
            obj.hessianEvalNo  = aux + 1;
        end
        
        function decrementHessianEvalNo(obj)
        %Decrements hessianEvalNo by one
            aux = obj.hessianEvalNo;
            obj.hessianEvalNo  = aux - 1;
        end
        
        function incrementBy(obj, FunctionGradientHessian)
        %Increments functionEvalNo by FunctionGradientHessian(1),
        %gradientEvalNo by FunctionGradientHessian(2)
        %hessianEvalNo by FunctionGradientHessian(3)
            auxF = obj.functionEvalNo;
            auxG = obj.gradientEvalNo;
            auxH = obj.hessianEvalNo;
            obj.functionEvalNo  = auxF + FunctionGradientHessian(1);
            obj.gradientEvalNo  = auxG + FunctionGradientHessian(2);
            obj.hessianEvalNo  = auxH + FunctionGradientHessian(3);
        end
        
        function str = DISPLAY(EN)
            str = strcat('functionEvalNo: ', num2str(EN.functionEvalNo), ' gradientEvalNo: ', num2str(EN.gradientEvalNo), ' hessianEvalNo: ', num2str(EN.hessianEvalNo));
        end
    end
    
end

