function [ outVal, outGr, outHes ] = ExtTridiag2( x0, VGH)
% Extended Tridiagonal 2 function
    
    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    c = 0.1;
    
    % computes the value of function in point x0
    if VGH(1) > 0
        for i = 1:n-1
            outVal = double(outVal + (x0(i)*x0(i+1)-1).^2 + c*(x0(i)+1)*(x0(i+1)+1));
        end
    end;
    
    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr(1) =  double(c*(1 + x0(2)) + 2*x0(2)*(-1 + x0(1)*x0(2)));
        for i=2:n-1
            outGr(i) = double(c*(1+x0(i-1))+2*x0(i-1)*(-1+x0(i-1)*x0(i))+c*(1+x0(i+1))+2*x0(i+1)*(-1+x0(i)*x0(i+1)));
        end
        outGr(n) = double(c*(1 + x0(n-1)) + 2*x0(n-1)*(-1 + x0(n-1)*x0(n)));
    end;
    
    % computes the numerical values of Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n,n);
        % computes first row
        outHes(1,1) =  double(2*x0(2).^2);
        outHes(1,2) = double(c+2*x0(1)*x0(2)+2*(-1+x0(1)*x0(2)));
        % computes rows from 2 to dim-1
        for i = 2:n-1
            outHes(i,i-1) = double(c+2*x0(i-1)*x0(i)+2*(-1+x0(i-1)*x0(i)));
            outHes(i,i) = double(2*x0(i-1).^2+2*x0(i+1).^2);
            outHes(i,i+1) = double(c+2*x0(i)*x0(i+1)+2*(-1+x0(i)*x0(i+1)));
        end
        % computes last row 
        outHes(n,n) = double(2*x0(n-1).^2);
        outHes(n,n-1) = double(c+2*x0(n-1)*x0(n)+2*(-1+x0(n-1)*x0(n)));
    else
        outHes = 0;
    end;
       
end