function [ outVal, outGr, outHes ] = LIARWHDNew( x, VGH )
% LIARWHD (CUTE) function

	n = length(x);
	outVal = 0;
	outGr = 0;
	
	% computes the value of function in point x
	if VGH(1) > 0
        for i = 1:n
			outVal = outVal + 4*(x(i)^2 - x(1))^2 + (x(i)-1)^2;
		end
	end
	
	% computes the numerical gradient value of function in point x
	if VGH(2) > 0
		outGr = zeros(n,1);
        auxVal = sum((x(2:end)).^2 - 1);
        outGr(1) = 16*x(1)^3 - 24*x(1)^2 + 10*x(1) - 2 - 8*auxVal;
        for i = 2:n 
            outGr(i) = 16*(x(i)^3 - x(1)*x(i)) + 2*(x(i)-1);
        end
	end

	% computes the numerical Hessian of function in point x
	if VGH(3) > 0
		outHes = zeros(n,n);
        outHes(1,1) = 48*x(1)^2 - 48*x(1) + 10 + 8*(n-1);
        for i = 2:n
            outHes(1,i) = -16*x(i);
            outHes(i,i) = 48*x(i)^2 - 16*x(1)+2;
            outHes(i,1) = -16*x(i);
        end
	else
		outHes = 0;
	end
	
end
