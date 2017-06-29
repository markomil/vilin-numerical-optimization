function [ outVal, outGr, outHes ] = GenTridiag1( x0, VGH)
% Generalized Tridiagonal 1 function:
    
    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    
    % computes the value of function in point x0
    if VGH(1) > 0
        for i = 1:n-1
            outVal = double(outVal + (x0(i)+x0(i+1)-3)^2 + (x0(i)-x0(i+1)+1)^4);
        end
    end;
    
    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr(1) = double(2*(x0(1)+x0(2)-3)+4*(x0(1)-x0(2)+1)^3);
        for i=2:n-1
            outGr(i) = double(2*(x0(i-1)+x0(i)-3)-4*(x0(i-1)-x0(i)+1)^3 + 2*(x0(i)+x0(i+1)-3)+4*(x0(i)-x0(i+1)+1)^3);
        end
        outGr(n) = double(2*(x0(n-1)+x0(n)-3)-4*(x0(n-1)-x0(n)+1)^3);
    end;
    
    % computes the numerical values of Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n,n);
        outHes(1,1) = 2 + 12*(x0(1)-x0(2)+1)^2;
        outHes(1,2) = 2 - 12*(x0(1)-x0(2)+1)^2;
        for i = 2:n-1
            outHes(i,i-1) = 2 - 12*(x0(i-1)-x0(i)+1)^2; 
            outHes(i,i) = 2 + 12*(x0(i-1)-x0(i)+1)^2 + 2 + 12*(x0(i)-x0(i+1)+1)^2;
            outHes(i,i+1) = 2 - 12*(x0(i)-x0(i+1)+1)^2;
        end
        outHes(n,n-1) = 2 - 12*(x0(n-1)-x0(n)+1)^2;
        outHes(n,n) = 2 + 12*(x0(n-1)-x0(n)+1)^2;
    else
        outHes = 0;
    end;
       
end
