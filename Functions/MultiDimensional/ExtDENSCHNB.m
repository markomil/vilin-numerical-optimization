function [outVal, outGr, outHes ] = ExtDENSCHNB( x, VGH )
% Extended DENSCHNB (CUTE) function

    n = length(x);    
    outVal = 0;
    outGr = 0;
    outHes = 0;
    
    % computes the value of function in point x0
    if VGH(1) > 0
                
        for i = 1:(n/2)
            outVal = outVal + (x(2*i-1)-2)^2+(x(2*i-1)-2)^2*x(2*i)^2+(x(2*i)+1)^2;
        end
    end
    
    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr = zeros(n,1);

        for i = 1:n
            if (mod(i,2)==1)
                outGr(i) = 2*(-2 + x(i))*(1 + x(i + 1)^2);
            else
                outGr(i) = 2*(1 + x(i) + (-2 + x(i-1))^2*x(i));
            end
        end
    end

    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n,n); 

        for i = 1:(n/2)   
            outHes(2*i-1,2*i-1) = 2*(1+x(2*i)^2);
            outHes(2*i,2*i) = 2*(1+(x(2*i-1)-2)^2);
            outHes(2*i,2*i-1) = 4*(-2 + x(2*i-1))*x(2*i);
            outHes(2*i-1,2*i) = outHes(2*i,2*i-1);
        end
    end
    
end

