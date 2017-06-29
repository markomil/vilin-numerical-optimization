function [ outVal, outGr, outHes ] = NONSCOMP( x0, VGH)
% NONSCOMP CUTE function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    c = 4;
    
    % computes the value of function in point x0
    if VGH(1) > 0
        outVal = (x0(1)-1)^2;
        for i = 2:n
            outVal = outVal + c*(x0(i)-(x0(i-1))^2)^2;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr(1) = 2*(-1 + x0(1)) - 4*c*x0(1)*(x0(2) - x0(1)^2);
        for i = 2:n-1
            outGr(i) = 2*c*(x0(i)-x0(i-1)^2) - 4*c*x0(i)*(x0(i+1) -x0(i)^2);
        end;
        outGr(n) = 2*c*(x0(n) - x0(n-1)^2);
    end;
     
    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        outHes(1, 1) = 2 + 8*c*x0(1)^2 - 4*c*(x0(2) - x0(1)^2);
        for i = 2:n
            outHes(i, i-1) = -4*c*x0(i-1);
            outHes(i-1, i) = -4*c*x0(i-1);
            if i == n
                outHes(i, i) = 2*c;
            else
                outHes(i, i) = 2*c + 8*c*x0(i)^2 - 4*c*(-x0(i)^2 + x0(i+1));
            end;
        end
    else
        outHes = 0;
    end;
  
end

