function [ outVal, outGr, outHes ] = PartPertQuad( x0, VGH )
% Partial Perturbed Quadratic function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	outHes = 0;
    
	%constants
	c = 100;

    if VGH(1) > 0
		outVal = x0(1)^2;
		innerSum = 0;
        
		for i = 1:n
			innerSum = innerSum + x0(i);
			outVal = outVal + i*x0(i)^2 + (1/c)*innerSum^2;
		end
    end
    
    if VGH(2) > 0
        outGr = zeros(n,1);
        auxVec = n:-1:1;
        auxSum = sum(auxVec.*x0);
        outGr(1) = 4*x0(1) + 2/c*auxSum;
        
        for i = 2:n
            auxSum = auxSum - sum(x0(1:i-1));
            outGr(i) = 2*i*x0(i) + 2/c*auxSum;
        end
    end
   
    if VGH(3) > 0
		outHes = zeros(n,n);
        
        for i = 1:n-1
            outHes(i, i) = 2/c * (n-i+1) + 2*i;
			for j = i+1:n
                outHes(i, j) = 2/c * (n-j+1);
				outHes(j, i) = outHes(i,j);
			end
        end
        outHes(1, 1) = 4 + 2/c*n;
        outHes(n, n) = 2/c + 2*n;
    end
    
end
