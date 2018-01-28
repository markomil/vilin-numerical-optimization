classdef LineSearchParams
    %Holds params for line search methods
    
    properties
        beta, rho, m, sigma, ksi, vals, grad, dir, startingPoint, tInitStart, tPrev, eps,
        it, theta, gamma, w, C
    end
    
    methods
        function init = LineSearchParams(methodParams, vals, grad, dir, startingPoint, tPrev, it, C)
             init.beta = methodParams.beta; 
             init.rho = methodParams.rho;
             init.m = methodParams.m;
             init.sigma = methodParams.sigma;
             init.ksi = methodParams.ksi;
             init.vals = vals;
             init.grad = grad;
             init.dir = dir;
             init.startingPoint = startingPoint;
             init.tInitStart = methodParams.startingPoint;
             init.tPrev = tPrev; %TODO add this in order to have prev line
             %search for two heuristic methods
             init.eps = methodParams.epsilon;
             init.it = it;
             init.theta = 0.5; %  used in the update rule in ApproxWolfe
             init.gamma = 0.66; % determines when a bisection step is performed in ApproxWolfe
             init.w = 1e-3; % used in switching from Wolfe to approximate Wolfe condition
             init.C = C;
             if nargin < 8
                 init.C = abs(vals(end));
             end
        end
    end
    
end

