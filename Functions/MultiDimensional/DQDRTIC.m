function [ outVal, outGr, outHes ] = DQDRTIC( x0, VGH )
% DQDRTIC (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	outHes = 0;

	%constants
	c = 100;
	d = 100;

	if VGH(1) > 0
		for i=1:n-2
			outVal = outVal + (x0(i)^2 + c*x0(i+1)^2 + d*x0(i+2)^2);
		end
	end
	
	if VGH(2) > 0
		outGr = zeros(n,1);
		for i=1:n-2
			outGr(i) = outGr(i) + 2*x0(i);
			outGr(i+1) = outGr(i+1) + 2*c*x0(i+1);
			outGr(i+2) = outGr(i+2) + 2*d*x0(i+2);
		end
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
		outHes(1,1) = 2;
		outHes(2,2) = 2 + 2*c;

		for i=3:n-2
			outHes(i,i) = 2 + 2*c + 2*d;
		end

		outHes(n-1, n-1) = 2*c + 2*d;
		outHes(n, n) = 2*d;
    end
        
end

