function [ outVal, outGr, outHes ] = FLETCHCR( x0, VGH )
% FLETCHCR (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	% constants
	c = 100;

	if VGH(1) > 0
		for i=1:n-1
			outVal = outVal + (x0(i+1)-x0(i) + 1 - x0(i)^2)^2;
		end
		outVal = c*outVal;
	end
	
	if VGH(2) > 0
		outGr = zeros(n,1);
		outGr(1) = 2*c*(-2*x0(1)-1)*(x0(2)-x0(1) + 1 - x0(1)^2);
		for i=2:n-1
			outGr(i) = 2*c*(x0(i)-x0(i-1) + 1 - x0(i-1)^2) + 2*c*(-2*x0(i)-1)*(x0(i+1)-x0(i) + 1 - x0(i)^2);
		end
		outGr(n) = 2*c*(x0(n)-x0(n-1) + 1 - x0(n-1)^2);
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
		outHes(1,1) = 2*c*(-1-2*x0(1))^2 - 4*c*(1-x0(1)-x0(1)^2+x0(2));
		outHes(1,2) = 2*c*(-1-2*x0(1));
		for i=2:n-1
			outHes(i,i-1) = 2*c*(-1-2*x0(i-1));
			outHes(i,i)   = 2*c - 4*c*(x0(i+1)-x0(i) + 1 - x0(i)^2) + 2*c*(-1-2*x0(i))^2;
			outHes(i,i+1) = 2*c*(-2*x0(i)-1);
		end
		outHes(n,n-1) = 2*c*(-1-2*x0(n-1));
		outHes(n,n)   = 2*c;
    else
        outHes = 0;
    end;
    
end

