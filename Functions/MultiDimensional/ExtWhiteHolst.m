function [ outVal, outGr, outHes ] = ExtWhiteHolst( x0, VGH )
% Extended White and Holst function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    c = 100;
    
    % computes the value of function in point x0
    if VGH(1) > 0
        i=1;
        while (i<=n/2)
            outVal = outVal + c*(x0(2*i)-(x0(2*i-1)).^3).^2 + (1-x0(2*i-1)).^2;
            i=i+1;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        if mod(n,2) == 1
            n = n-1;
        end;
        for i = 1:n
            if mod(i,2) == 0
                outGr(i) = c*(2*x0(i)-2*(x0(i-1)).^3);
            else
                outGr(i) = c*(6*(x0(i)).^5 -6*x0(i+1)*(x0(i)).^2)+(2*x0(i)-2);
            end;
        end;
    end;
    
	% computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            if mod(i,2) == 0
                outHes(i, i-1) = double(-6*c*x0(i-1)^2);
                outHes(i, i) = 2*c;
            else
                outHes(i, i) = double(18*c*x0(i)^4 - 12*c*x0(i)*(x0(i+1) - x0(i)^3) + 2); 
                outHes(i, i+1) = double(-6*c*x0(i)^2);
            end;
        end;
    else
        outHes = 0;
    end;
    
end

