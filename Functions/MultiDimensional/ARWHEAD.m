function [ outVal, outGr, outHes ] = ARWHEAD( x, VGH )
% ARWHEAD (CUTE) function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    c = 0.5;
    
    if VGH(1) > 0
        for i = 1:n-1
            outVal = outVal + -4*x(i)+3 + (x(i)^2 + x(n)^2)^2;
        end
    end

    if VGH(2) > 0
        inner = 0;
        for i = 1:n-1
            outGr(i) = 4*(x(i)^3 + x(i)*x(n)^2 - 1);
            inner = inner + x(i)^2;
        end
        outGr(n) = 4*x(n)*(inner + 2*x(n)^2);
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        inner = 0;
        for i = 1:n-1
            outHes(i, i) = 4 * (3*x(i)^2 + x(n)^2);
            outHes(n, i) = 8 * x(i) * x(n);
            outHes(i, n) = 8 * x(i) * x(n);
            inner = inner + x(i)^2;
        end
        outHes(n, n) = 4 * (inner + 6*x(n)^2);
    else
        outHes = 0;
    end

end