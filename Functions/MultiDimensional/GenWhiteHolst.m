function [ outVal, outGr, outHes ] = GenWhiteHolst( x0, VGH )
% Generalized White and Holst function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    c = 100;
    
    % computes the value of function in point x0
    if VGH(1) > 0
        i=1;
        while (i<=n-1)
            outVal = outVal + c*(x0(i+1)-(x0(i))^3)^2 + (1-x0(i))^2;
            i=i+1;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr(1) = -2*(1-x0(1)) - 6*c*x0(1)^2*(-x0(1)^3 + x0(2));
        for i = 2:n-1
            outGr(i) = -2*(1-x0(i)) - 6*c*x0(i)^2*(-x0(i)^3 + x0(i+1)) + 2*c*(-x0(i-1)^3 + x0(i));
        end;
        outGr(n) = 2*c*(-x0(n-1)^3 + x0(n));
    end;
    
	% computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n-1
            outHes(i, i+1) = double(-6*c*x0(i)^2);
            outHes(i+1, i) = double(-6*c*x0(i)^2);
            outHes(i, i) = double(2 + 2*c + 18*c*x0(i)^4 - 12*c*x0(i)*(-x0(i)^3 + x0(i+1)));
        end;
        outHes(1, 1) = outHes(1, 1) - 2*c;
        outHes(n, n) = 2*c;
    else
        outHes = 0;
    end;
    
end

