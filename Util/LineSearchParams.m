classdef LineSearchParams
    %Holds params for line search methods
    
    properties
        beta, rho, m, sigma, ksi, val, grad, dir, startingPoint, tStart, tInitStart, eps, theta, gama 
    end
    
    methods
        function init = LineSearchParams(methodParams, val, grad, dir, startingPoint, tStart)
             init.beta = methodParams.beta; 
             init.rho = methodParams.rho;
             init.m = methodParams.m;
             init.sigma = methodParams.sigma;
             init.ksi = methodParams.ksi;
             init.val = val;
             init.grad = grad;
             init.dir = dir;
             init.startingPoint = startingPoint;
             init.tStart = tStart;
             init.tInitStart = methodParams.startingPoint;
             init.eps = methodParams.epsilon;
             init.theta = 0.5;
             init.gama = 0.66;
         end
    end
    
end

