function [ outVal, outGr, outHes ] = QuadQF1( x, VGH )
% Quadratic QF1 function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + i*x(i)^2;
        end
        outVal = 0.5*outVal - x(n);
    end

    if VGH(2) > 0
        for i = 1:n
            outGr(i) = i*x(i);
        end
        outGr(n) = outGr(n) -1;
    end

    if VGH(3) > 0
        diagEl = 1:n;
        outHes = diag(diagEl);
    else
        outHes = 0;
    end

end
