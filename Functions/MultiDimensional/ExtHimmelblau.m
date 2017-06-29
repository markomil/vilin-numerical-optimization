function [outVal, outGr, outHes ] = ExtHimmelblau( x0, VGH)
% Extended Himmelblau function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
        
    % computes the value of function in point x0
    if VGH(1) > 0
        i = 1;
        while (i<=n/2)
            outVal = outVal + (x0(2*i-1)^2 + x0(2*i) - 11).^2 + (x0(2*i-1) + x0(2*i)^2 - 7).^2;
            i = i + 1;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        if mod(n,2) == 1
            n = n-1;
        end;
        for i = 1:n/2
            outGr(2*i-1) = 4*x0(2*i-1)^3 + 2*x0(2*i)^2 + 4*x0(2*i-1)*x0(2*i) - 42*x0(2*i-1) - 14;
            outGr(2*i) = 4*x0(2*i)^3 + 2*x0(2*i-1)^2 + 4*x0(2*i-1)*x0(2*i) - 26*x0(2*i) - 22;
        end
    end;
    
    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        if mod(n,2) == 1
            n = n-1;
        end;
        for i = 1:n/2
            outHes(2*i-1, 2*i-1) = 12*x0(2*i-1)^2 + 4*x0(2*i) - 42;
            outHes(2*i, 2*i-1) = 4*x0(2*i) + 4*x0(2*i-1);
            outHes(2*i-1, 2*i) = 4*x0(2*i) + 4*x0(2*i-1);
            outHes(2*i, 2*i) = 12*x0(2*i)^2 + 4*x0(2*i-1) - 22;
        end;
    else
        outHes = 0;
    end;

end

