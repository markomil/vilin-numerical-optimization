function [ outVal, outGr, outHes ] = HIMMELHNew( x, VGH )
% HIMMELH (CUTE) function 

	n = length(x);
	outVal = 0;
	outGr = 0;
	
    % computes the value of function in point x
    if VGH(1) > 0
		for i = 1:n/2
			outVal = outVal - 3*x(2*i-1) - 2*x(2*i) + 2 + x(2*i-1)^3 + x(2*i)^2;
		end
    end
	
    % computes the numerical gradient value of function in point x
    if VGH(2) > 0
		outGr = zeros(n,1);
        for i = 1:n
            if mod(i,2) == 0
                outGr(i) = -2 + 2*x(i);
            else
                outGr(i) = -3 + 3*x(i)^2;
            end
        end
    end

    % computes the numerical Hessian of function in point x
	if VGH(3) > 0
		outHes = zeros(n,n);
        for i = 1:n
            if mod(i,2) == 0
                outHes(i,i) = 2;
            else
                outHes(i,i) = 6*x(i);
            end
        end
	else
		outHes = 0;
	end
end
