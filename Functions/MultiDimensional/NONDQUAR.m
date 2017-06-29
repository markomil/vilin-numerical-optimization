function [ outVal, outGr, outHes ] = NONDQUAR( x0, VGH )
% NONDQUAR (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	if VGH(1) > 0
		outVal = (x0(1)-x0(2))^2 + (x0(n-1)+x0(n))^2;
		for i=1:n-2
			outVal = outVal + (x0(i)+x0(i+1)+x0(n))^4;
		end
	end
	
	if VGH(2) > 0
		outGr = zeros(n,1);
		outGr(1) = 2*(x0(1)-x0(2));
		outGr(2) = -2*(x0(1)-x0(2));
		for i=1:n-2
			ci = 4*(x0(i)+x0(i+1)+x0(n))^3;
			outGr(i) = outGr(i) + ci;
			outGr(i+1) = outGr(i+1) + ci;
			outGr(n) = outGr(n) + ci;
		end
		cn1 = 2*(x0(n-1)+x0(n));
		outGr(n-1) = 4*(x0(n-2)+x0(n-1)+x0(n))^3 + cn1;
		outGr(n) = outGr(n) + cn1;
    end

    if VGH(3) > 0
		outHes = zeros(n,n);

		outHes(1,1) = 2;
		outHes(1,2) = -2;
		outHes(2,1) = -2;
		outHes(2,2) = 2;

		for i=1:n-2
			ci = 12*(x0(i)+x0(i+1)+x0(n))^2;

			outHes(i,i) = outHes(i,i) + ci;
			outHes(i,i+1) = outHes(i,i+1) + ci;
			outHes(i,n) = outHes(i,n) + ci;

			outHes(i+1,i) = outHes(i+1,i) + ci;
			outHes(i+1,i+1) = outHes(i+1,i+1) + ci;
			outHes(i+1,n) = outHes(i+1,n) + ci;

			outHes(n,i) = outHes(n,i) + ci;
			outHes(n,i+1) = outHes(n,i+1) + ci;
			outHes(n,n) = outHes(n,n) + ci;
		end
		
		outHes(n-1,n-1) = outHes(n-1,n-1) + 2; 
		outHes(n-1,n) = outHes(n-1,n) + 2; 

		outHes(n,n-1) = outHes(n,n-1) + 2;
		outHes(n,n) = outHes(n,n) + 2;
    else
        outHes = 0;
    end
    
end

