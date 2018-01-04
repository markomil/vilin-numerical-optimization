function [ outVal, outGr, outHes ] = ExtHiebert( x0, VGH )
% Extended Hiebert function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	% constants
	c1 = 10;
	c2 = 50000;

    if VGH(1) > 0
		for i=1:n/2
			outVal = outVal + (x0(2*i-1)-c1)^2 + (x0(2*i-1)*x0(2*i)-c2)^2;
		end
    end
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		for i=1:n/2
			outGr(2*i-1) = 2*(x0(2*i-1)-c1) + 2*x0(2*i)*(x0(2*i-1)*x0(2*i)-c2);
			outGr(2*i) = 2*x0(2*i-1)*(x0(2*i-1)*x0(2*i)-c2);
		end
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
		for i=1:n/2
			outHes(2*i-1, 2*i-1) = 2 + 2*x0(2*i)^2;
			outHes(2*i-1, 2*i) = 4*x0(2*i-1)*x0(2*i)-2*c2;

			outHes(2*i, 2*i-1) = 4*x0(2*i-1)*x0(2*i)-2*c2;
			outHes(2*i, 2*i) = 2*x0(2*i-1)^2;
		end
	else
        outHes = 0;
    end;

end

