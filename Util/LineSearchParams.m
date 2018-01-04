classdef LineSearchParams
    %Holds params for line search methods
    
    properties
        beta, rho, m, sigma, ksi, vals, grad, dir, startingPoint, tInitStart, tPrev, eps, theta, gamma, it 
    end
    
    methods
        function init = LineSearchParams(methodParams, vals, grad, dir, startingPoint, tPrev, it)
             init.beta = methodParams.beta; 
             init.rho = methodParams.rho;
             init.m = methodParams.m;
             init.sigma = methodParams.sigma;
             init.ksi = methodParams.ksi;
             init.vals = vals;
             init.grad = grad;
             init.dir = dir;
             init.startingPoint = startingPoint;
             init.tPrev = tPrev; %TODO add this in order to have prev line
             %search for two heuristic methods
             init.tInitStart = methodParams.startingPoint;
             init.eps = methodParams.epsilon;
             init.theta = 0.5;
             init.gamma = 0.66;
             init.it = it;
         end
    end
    
end

