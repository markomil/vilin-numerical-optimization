classdef MathodParams

    properties
        starting_point, step_size, variables_no, max_iteration_no, epsilon, workPrec,
        step_size_min, od, do, beta, sigma, rho, m, ksi, nu, startingPoint, 
        lineSearchMethod
    end
    
    methods
        function obj = MathodParams(sp, ss, vn, min,  eps, workPrec, ssm, od, do, beta, sigma, rho, m, ksi, startingPoint, lineSearchMethod) 
           if nargin > 0
               obj.starting_point = sp;
               obj.step_size = ss;
               obj.variables_no = vn;
               obj.max_iteration_no = min;
               obj.epsilon = eps;
               obj.workPrec = workPrec;
               obj.step_size_min = ssm;
               obj.od = od;
               obj.do = do;
               obj.beta = beta;
               obj.sigma = sigma;
               obj.rho = rho;
               obj.ksi = ksi;
               obj.nu = 0.1; % threshold for restarting beta in CG methods
               obj.m = m;
               obj.startingPoint = startingPoint;
               obj.lineSearchMethod = lineSearchMethod;
           end
        end
    end
    
end

