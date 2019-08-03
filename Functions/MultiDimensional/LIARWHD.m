function [ outVal, outGr, outHes ] = LIARWHD( x, VGH )
% LIARWHD (CUTE) function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + 4*(-x(1)+x(i)^2)^2 + (x(i)-1)^2;
        end
    end

    if VGH(2) > 0
        c = 10 + (n-1)*8;
        outGr(1) = 16*x(1)^3 - 24*x(1)^2 + c*x(1) - 2;
        for i = 2:n
            outGr(1) = outGr(1) - 8*x(i)^2;
            outGr(i) = 2*(-8*x(1)*x(i) + 8*x(i)^3 + x(i) - 1);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        c = 10 + (n-1)*8;
        for i = 2:n
            outHes(i, i) = -16*x(1) + 48*x(i)^2 + 2;
            outHes(1, i) = -16*x(i);
            outHes(i, 1) = -16*x(i);
        end
        outHes(1, 1) = 48*x(1)^2 - 48*x(1) + c;
    else outHes = 0;
    end

end