function [ outVal, outGr, outHes ] = NewFunctionTemplate( x, VGH )
% ================================================================
% Template for adding new multidimensional function to Vilin.
% To add new function modify this file and save to 'Functions/MultiDimensional/'.
% To add starting point for new function check 'Util/StartingPointGenerator.m'.
% ================================================================

    n = length(x);
    outVal = 0;
    outGr = 0;
    outHes = 0;
    
    if VGH(1) == 1
        % compute function value, replace line below with actual code
        throw (MException('NumOpt:implementationError', 'Computing value for %s not implemented.', mfilename));
    end

    if VGH(2) == 1
        outGr = zeros(n, 1);
        % compute gradient, replace line below with actual code
        throw (MException('NumOpt:implementationError', 'Computing gradient for %s not implemented.', mfilename));
    end

    if VGH(3) == 1
        outHes = zeros(n, n);
        % compute hessian, replace line below with actual code
        throw (MException('NumOpt:implementationError', 'Computing hessian for %s not implemented.', mfilename));
    end

end
