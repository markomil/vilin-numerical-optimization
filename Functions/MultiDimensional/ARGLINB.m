function [ outVal, outGr, outHes ] = ARGLINB( x0, VGH )
% ARGLINB (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	innerSums = 0;

    if VGH(1) > 0
		innerSums = zeros(n,1);
		for i=1:n
			innerSum = 0;
			for j=1:n
				innerSum = innerSum + i*j*x0(j) - 1;
			end
			outVal = outVal + innerSum^2;
			innerSums(i) = innerSum;
		end
    end
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		for i=1:n
			innerSum = getInnerSum(i, innerSums, n, x0, VGH);	
			for j=1:n
				outGr(j) = outGr(j) + 2*i*j*innerSum;
			end
		end
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
        for i=1:n
            for j=1:n
				for k=1:n
					outHes(j, k) = outHes(j, k) + 2*i^2*k*j;
				end
            end
        end
    else
        outHes = 0;
    end;
    
end
            
function innerSum = getInnerSum(i, innerSums, n, x0, VGH)
    
    if VGH(1) > 0 %innerSums is initialized
        innerSum = innerSums(i);
    else
        innerSum = 0;
        for j=1:n
            innerSum = innerSum + i*j*x0(j) - 1;
        end
    end
    
end

