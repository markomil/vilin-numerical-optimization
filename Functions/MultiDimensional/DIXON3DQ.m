function [outVal, outGr, outHes ] = DIXON3DQ( x0, VGH)
% DIXON3DQ (CUTE) function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
        
    % computes the value of function in point x0
    if VGH(1) > 0
        outVal = (x0(1) - 1)^2 + (x0(n) - 1)^2;
        for i = 1:n-1
            outVal = outVal + (x0(i) - x0(i+1))^2;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr(1) = 2*(-1 + x0(1)) + 2*(x0(1) - x0(2));
        for i = 2:n-1
            outGr(i) = -2*(x0(i-1) - x0(i)) + 2*(x0(i) - x0(i+1));
        end;
        outGr(n) = -2*(x0(n-1) - x0(n)) + 2*(-1 + x0(n));
    end;
       
    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        outHes(1, 1) = 4;
        for i = 2:n
            outHes(i, i-1) = -2;
            outHes(i-1, i) = -2;
            outHes(i, i) = 4;
        end
    else
        outHes = 0;
    end;

end

