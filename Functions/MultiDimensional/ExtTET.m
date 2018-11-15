function [ outVal,outGr,outHes ] = ExtTET( x0,VGH )
% Extended TET (Three exponential terms) function

    n = length(x0);
    assert (mod(n,2)==0)
    c = 0.90483741803596; %exp(-0.1)

    outVal = 0;
    outGr = zeros(n, 1);

   	% computes the value of function in point x0
    if VGH(1) > 0
       for i=2:2:n
           outVal = outVal + exp(x0(i-1)+3*x0(i)-0.1)+exp(x0(i-1)-3*x0(i)-0.1)+exp(-x0(i-1)-0.1);
       end
    end

    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
       for i=1:n
          if mod(i,2) == 1
              outGr(i) = c*((exp(x0(i) - 3*x0(i+1)) + exp(x0(i) + 3*x0(i+1)))*exp(x0(i)) - 1)*exp(-x0(i));
          else
              outGr(i) = -3*c*exp(x0(i-1) - 3*x0(i)) + 3*c*exp(x0(i-1) + 3*x0(i));
          end
       end
    end

    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
       outHes = zeros(n, n);     
       for i=2:2:n
           outHes(i-1,i-1) = c*((exp(x0(i-1) - 3*x0(i)) + exp(x0(i-1) + 3*x0(i)))*exp(x0(i-1)) + 1)*exp(-x0(i-1));
           outHes(i,i)= 9*c*exp(x0(i-1) - 3*x0(i)) + 9*c*exp(x0(i-1) + 3*x0(i));
           outHes(i,i-1) = -3*c*exp(x0(i-1) - 3*x0(i)) + 3*c*exp(x0(i-1) + 3*x0(i));
           outHes(i-1,i) = -3*c*exp(x0(i-1) - 3*x0(i)) + 3*c*exp(x0(i-1) + 3*x0(i));
       end
    else
        outHes = 0;
    end
    
end

