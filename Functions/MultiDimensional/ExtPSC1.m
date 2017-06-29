function [ outVal, outGr, outHes ] = ExtPSC1( x0, VGH )
% Extended PSC1 function

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    
    
    % computes the value of function in point x0
    if VGH(1) > 0
        i=1;
        while (i<=n/2)
            outVal = outVal + (x0(2*i-1)^2 + x0(2*i)^2 + x0(2*i-1)*x0(2*i))^2 + ...
                sin(x0(2*i-1))^2 + cos(x0(2*i))^2; 
            i=i+1;
        end;
    end;

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        for i = 1:n/2
            outGr(2*i-1) = 2*(2*x0(2*i-1)+x0(2*i))*(x0(2*i-1)^2+x0(2*i-1)*x0(2*i)+x0(2*i)^2) + ...
                    2*cos(x0(2*i-1))*sin(x0(2*i-1));
            outGr(2*i) = 2*(x0(2*i-1)+2*x0(2*i))*(x0(2*i-1)^2+x0(2*i-1)*x0(2*i)+x0(2*i)^2) - ...
                    2*cos(x0(2*i))*sin(x0(2*i));
        end;
    end;
    
     
  	% computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(n, n);
        for j = 1:n/2
            i = 2*j-1;
            % diagonal elements
            outHes(i, i) = double(2*(2*x0(i)+x0(i+1))^2 + 4*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2)) + ...
                        2*cos(x0(i))^2 - 2*sin(x0(i))^2;
            outHes(i+1, i+1) = double(2*(x0(i)+2*x0(i+1))^2 + 4*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2)) - ...
                2*cos(x0(i+1))^2 + 2*sin(x0(i+1))^2;
            
            % non diagonal elements
            outHes(i, i+1) = double(2*(2*x0(i)+x0(i+1))*(x0(i)+2*x0(i+1)) + 2*(x0(i)^2+x0(i)*x0(i+1)+x0(i+1)^2));
            outHes(i+1, i) = outHes(i, i+1);
                   
        end;
    else
        outHes = 0;
    end;
  
end

