function [ outVal, outGr, outHes ] = EDENSCH( x, VGH )
% EDENSCH (CUTE) function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        outVal = 16;
        for i = 1:n-1
            outVal = outVal + (x(i)-2)^4 + (x(i)*x(i+1)-2*x(i+1))^2 + (x(i+1)+1)^2;
        end
    end

    if VGH(2) > 0
        outGr(1) = 2*(x(1)-2)*(2*x(1)^2 - 8*x(1) + x(2)^2 + 8);
        outGr(n) = 2*(x(n-1)^2 - 4*x(n-1) + 5)*x(n) + 2;
        for i = 2:n-1
            outGr(i) = 2*( ((x(i-1)-2)^2)*x(i) + (x(i)-2)*x(i+1)^2 + 2*(x(i)-2)^3 + x(i) + 1);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        outHes(1, 1) = 2* (6*x(1)^2 - 24*x(1) + x(2)^2 + 24);
        outHes(1, 2) = 4*(x(1)-2)*x(2);
        outHes(n, n-1) = 4*(x(n-1)-2)*x(n);
        outHes(n, n) = 2* (x(n-1)^2 - 4*x(n-1) + 5);
        for i = 2:n-1
            outHes(i, i-1) = 4*(x(i-1)-2)*x(i);
            outHes(i, i) = 2*( (x(i-1)-2)^2 + 6*(x(i)-2)^2 + x(i+1)^2 + 1);
            outHes(i, i+1) = 4*(x(i)-2)*x(i+1);
        end
    else outHes = 0;
    end

end