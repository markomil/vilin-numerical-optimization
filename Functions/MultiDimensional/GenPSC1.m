function [ outVal, outGr, outHes ] = GenPSC1( x0, VGH )
% Generalized PSC1 function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
        
    % computes the value of function in point x0
    if VGH(1) > 0
        i=1;
        while (i<=n-1)
            outVal = outVal + (x0(i)^2 + x0(i+1)^2 + x0(i)*x0(i+1))^2 + sin(x0(i))^2 + cos(x0(i+1))^2; 
            i=i+1;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr(1) = 2*(2*x0(1)+x0(2))*(x0(1)^2+x0(1)*x0(2)+x0(2)^2)+2*cos(x0(1))*sin(x0(1));
        for i = 2:n-1
            outGr(i) = 2*(x0(i-1)+2*x0(i))*(x0(i-1)^2+x0(i-1)*x0(i)+x0(i)^2) + ...
                       2*(2*x0(i)+x0(i+1))*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2);
        end;
        outGr(n) = 2*(x0(n-1)+2*x0(n))*(x0(n-1)^2+x0(n-1)*x0(n)+x0(n)^2)-2*cos(x0(n))*sin(x0(n));
    end;
    
  	% computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n-1
            outHes(i, i+1) = double(2*(2*x0(i)+x0(i+1))*(x0(i)+2*x0(i+1)) + 2*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2));
            outHes(i+1, i) = double(2*(2*x0(i)+x0(i+1))*(x0(i)+2*x0(i+1)) + 2*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2));
            if i == 1
                outHes(i, i) = double(2*(2*x0(1)+x0(2))^2 + 4*(x0(1)^2+x0(1)*x0(2)+x0(2)^2)) + 2*cos(x0(1))^2 - 2*sin(x0(1))^2;
            else
                outHes(i, i) = double(2*(x0(i-1)+2*x0(i))^2 + 4*(x0(i-1)^2+x0(i-1)*x0(i)+x0(i)^2) + ... 
                    2*(2*x0(i)+x0(i+1))^2 + 4*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2));
            end;
        end;
        outHes(n, n) = double(2*(x0(n-1)+2*x0(n))^2 + 4*(x0(n-1)^2+x0(n-1)*x0(n)+x0(n)^2)) - 2*cos(x0(n))^2 + 2*sin(x0(n))^2;
    else
        outHes = 0;
    end;
    
end

