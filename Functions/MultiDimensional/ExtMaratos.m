function [ outVal, outGr, outHes ] = ExtMaratos( x0, VGH )
% Extended Maratos function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    c = 100;
    
    
    % computes the value of function in point x0
    if VGH(1) > 0
        i=1;
        while (i<=n/2)
            outVal = outVal + x0(2*i-1) + c*(x0(2*i-1)^2 + x0(2*i)^2-1)^2;
            i=i+1;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        for i = 1:n/2
            j = 2*i-1;
            outGr(j) = 1 + 4*c*x0(j)*(-1 + x0(j)^2 + x0(j+1)^2);
            outGr(j+1) = 4*c*x0(j+1)*(-1 + x0(j)^2 + x0(j+1)^2);
        end;
    end;
   
  	% computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        for j = 1:n/2
            i = 2*j-1;
            % diagonal elements
            outHes(i, i) = double(8*c*x0(i)^2 + 4*c*(-1+x0(i)^2+x0(i+1)^2));
            outHes(i+1, i+1) = double(8*c*x0(i+1)^2 + 4*c*(-1+x0(i)^2+x0(i+1)^2));
           
            % non diagonal elements
            outHes(i, i+1) = double(8*c*x0(i)*x0(i+1));
            outHes(i+1, i) = outHes(i, i+1);
                   
        end;
    else
        outHes = 0;
    end;
  
end

