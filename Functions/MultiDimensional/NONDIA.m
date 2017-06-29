function [ outVal, outGr, outHes ] = NONDIA( x0, VGH )
% NONDIA (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	%constants
	c = 100;

    if VGH(1) > 0
		outVal = (x0(1)-1)^2;
		for i = 2:n
			outVal = outVal + c*(x0(1)-x0(i-1)^2)^2;
		end
    end
    
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		outGr(1) = 2*(x0(1)-1) + 2*c*(1-2*x0(1))*(x0(1)-x0(1)^2); % + 2*c*(x0(1)-x0(2)^2);
		for i = 3:n
			outGr(1) = outGr(1) + 2*c*(x0(1)-x0(i-1)^2);
			outGr(i-1) = -4*c*x0(i-1)*(x0(1)-x0(i-1)^2);
		end
    end
        
    if VGH(3) > 0
        outHes = zeros(n,n);
        outHes(1,1) = 2 - 12*c*x0(1) + 12*c*x0(1)^2 + 2*(n-1)*c;
        for i = 3:n
			outHes(1,i-1) = outHes(1,i-1) - 4*c*x0(i-1);
			outHes(i-1,1) = outHes(i-1,1) - 4*c*x0(i-1);
			outHes(i-1,i-1) = outHes(i-1,i-1) + 8*c*x0(i-1)^2 - 4*c*(x0(1)-x0(i-1)^2);
        end
    else
        outHes = 0;
    end
 
end

