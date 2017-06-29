function [ outVal, outGr, outHes ] = AlmostPertQuad( x0, VGH )
% Almost Perturbed Quadratic function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	%constants
	c = 100;

    if VGH(1) > 0
		innerPart = (1/c)*(x0(1)+x0(n))^2;
		for i=1:n
			outVal = outVal + i*x0(i)^2 + innerPart;
		end
    end
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		innerPartDer = (2/c)*(x0(1)+x0(n));%derivative is the same for both vars
		for i=1:n
			outGr(i) = outGr(i) + 2*i*x0(i);
		end
		outGr(1) = outGr(1) + n*innerPartDer;
		outGr(n) = outGr(n) + n*innerPartDer;
    end

    if VGH(3) > 0
		outHes = zeros(n,n);

		for i=1:n
			outHes(i,i) = outHes(i,i) + 2*i;
		end
		outHes(1,1) = outHes(1,1) + (2*n)/c;
		outHes(1,n) = outHes(1,n) + (2*n)/c;
			
		outHes(n,1) = outHes(n,1) + (2*n)/c;
		outHes(n,n) = outHes(n,n) + (2*n)/c;
    else
        outHes = 0;
    end;
    
end
