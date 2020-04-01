function [ outVal, outGr, outHes ] = ENGVAL1( x, VGH )
% ENGVAL1 (CUTE) function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n-1
            outVal = outVal + (x(i)^2 + x(i+1)^2)^2 + (-4*x(i) + 3);
        end
    end

    if VGH(2) > 0
        outGr(1) = 4 * (x(1)^3 + x(1)*x(2)^2 - 1);
        outGr(n) = 4*x(n) * (x(n-1)^2 + x(n)^2);
        for i = 2:n-1
            outGr(i) = 4 * (x(i)*x(i-1)^2 + 2*x(i)^3 + x(i)*x(i+1)^2 - 1);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        outHes(1, 1) = 4 * (3*x(1)^2 + x(2)^2);
        outHes(1, 2) = 8 * x(1) * x(2);
        outHes(n, n-1) = 8 * x(n-1) * x(n);
        outHes(n, n) = 4 * (3*x(n)^2 + x(n-1)^2);
        for i = 2:n-1
            outHes(i, i-1) = 8*x(i-1)*x(i);
            outHes(i, i) = 4*(x(i-1)^2 + 6*x(i)^2 + x(i+1)^2);
            outHes(i, i+1) = 8*x(i)*x(i+1);
        end
    else
        outHes = 0;
    end

end