function [ outVal, outGr, outHes ] = EG2( x0, VGH )
% EG2 (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	outHes = 0;

	if VGH(1) > 0
		for i = 1:n-1
			outVal = outVal + sin(x0(1) + x0(i)^2 - 1) + 0.5*sin(x0(n)^2);
		end
	end
	
	if VGH(2) > 0
		outGr = zeros(n,1);
		outGr(1) = (2*x0(1)+1)*cos(x0(1) + x0(1)^2 - 1);
		outGr(n) = x0(n)*cos(x0(n)^2);

		for i = 2:n-1
			outGr(1) = outGr(1) + cos(x0(1) + x0(i)^2 - 1); 
			outGr(i) = outGr(i) + 2*x0(i)*cos(x0(1) + x0(i)^2 - 1);
			outGr(n) = outGr(n) + x0(n)*cos(x0(n)^2);
		end
	end

	if VGH(3) > 0
		outHes = zeros(n,n);
		outHes(1,1) = 2*cos(x0(1) + x0(1)^2 - 1) - (2*x0(1)+1)^2*sin(x0(1) + x0(1)^2 - 1);
		outHes(n,n) = cos(x0(n)^2) - 2*x0(n)^2*sin(x0(n)^2);

		for i = 2:n-1
			outHes(1,1) =  outHes(1,1) - sin(x0(1) + x0(i)^2 - 1);
			outHes(1,i) =  - 2*x0(i)*sin(x0(1) + x0(i)^2 - 1);

			outHes(i,1) = - 2*x0(i)*sin(x0(1) + x0(i)^2 - 1);
			outHes(i,i) = 2*cos(x0(1) + x0(i)^2 - 1) - 4*x0(i)^2*sin(x0(1) + x0(i)^2 - 1);
			
			outHes(n,n) = outHes(n,n) + cos(x0(n)^2) - 2*x0(n)^2*sin(x0(n)^2);
		end
	end
end
