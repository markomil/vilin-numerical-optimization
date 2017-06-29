function [ outVal, outGr, outHes ] = QUARTC( x0, VGH)
% QUARTC CUTE function 

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
        
    % computes the value of function in point x0
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + (x0(i)-1)^4;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        for i = 1:n
            outGr(i) = 4*(x0(i) - 1)^3;
        end;
    end;
     
    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            outHes(i, i) = 12*(x0(i) - 1)^2;
        end;
    else
        outHes = 0;
    end;
  
end



