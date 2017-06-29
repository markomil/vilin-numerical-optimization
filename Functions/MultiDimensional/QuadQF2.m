function [ outVal, outGr, outHes ] = QuadQF2( x0, VGH )
% Quadratic QF2 function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	% constants
	c = 0.5;

    if VGH(1) > 0
		for i = 1:n
			outVal = outVal + i*(x0(i)^2-1)^2;
		end
		outVal = c*outVal - x0(n);
    end
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		for i = 1:n
			outGr(i) = 4*c*i*x0(i)*(x0(i)^2-1);
		end
		outGr(n) = outGr(n) - 1;
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
		for i = 1:n
			outHes(i, i) = 12*c*i*x0(i)^2 - 4*c*i;
		end
	else
        outHes = 0;
    end
    
end

