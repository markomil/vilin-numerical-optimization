function [ outVal, outGr, outHes ] = Diagonal9( x, VGH )
% Diagonal 9 function

	n = length(x);
	outVal = 0;
	outGr = zeros(n, 1);
	c = 10000;
	
	% computes the value of function in point x
	if VGH(1) > 0
        for i = 1:n-1
			outVal = outVal + exp(x(i)) - i*x(i);
        end
        outVal = outVal + c*x(n)^2;
	end
	
	% computes the numerical gradient value of function in point x
	if VGH(2) > 0
		for i = 1:n-1
            outGr(i) = exp(x(i)) - i;
        end
        outGr(n) = 2*c*x(n);
	end

	% computes the numerical Hessian of function in point x0
    if VGH(3) > 0
		outHes = zeros(n,n);
        for i = 1:n-1
            outHes(i,i) = exp(x(i));
        end
        outHes(n,n) = 2*c;
    else
        outHes = 0;
    end;
    
end
