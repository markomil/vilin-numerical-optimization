function [ outVal, outGr, outHes ] = ExtTridiag1( x0, VGH)
% Extended Tridiagonal 1 function:
    
    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    
    % computes the value of function in point x0
    if VGH(1) > 0
        for i = 1:n/2
            outVal = outVal + (x0(2*i-1)+x0(2*i)-3)^2 + (x0(2*i-1)-x0(2*i)+1)^4;
        end
    end;
    
    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        for i=1:n/2
            outGr(2*i-1) = 2*(x0(2*i-1)+x0(2*i)-3)+4*(x0(2*i-1)-x0(2*i)+1)^3;
            outGr(2*i) = 2*(x0(2*i-1)+x0(2*i)-3)-4*(x0(2*i-1)-x0(2*i)+1)^3;
        end
    end;
    
    % computes the numerical values of Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n,n);
        for i = 1:n/2
            outHes(2*i-1,2*i-1) = 2 + 12*(x0(2*i-1)-x0(2*i)+1)^2;
            outHes(2*i-1,2*i) = 2 - 12*(x0(2*i-1)-x0(2*i)+1)^2;
            outHes(2*i,2*i-1) = 2 - 12*(x0(2*i-1)-x0(2*i)+1)^2;
            outHes(2*i,2*i) = 2 + 12*(x0(2*i-1)-x0(2*i)+1)^2;
        end
    else
        outHes = 0;
    end;
       
end