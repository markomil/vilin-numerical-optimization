function [ outVal, outGr, outHes ] = POWER( x, VGH )
% POWER (CUTE) function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + (i * x(i))^2;
        end
    end

    if VGH(2) > 0
        for i = 1:n
            outGr(i) = 2*i^2 * x(i);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            outHes(i, i) = 2*i^2;
        end
    else outHes = 0;
    end

end