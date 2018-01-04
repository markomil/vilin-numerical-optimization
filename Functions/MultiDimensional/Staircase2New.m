function [ outVal, outGr, outHes ] = Staircase2New( x0, VGH )
% Staircase 2 function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	innerSums = zeros(n,1);

    if VGH(1) > 0
		for i = 1:n
			innerSum = -i;
            innerSum = innerSum + sum(x0(1:i));
			
			innerSums(i) = innerSum;
			outVal = outVal + innerSum^2;
		end
    end
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		for i = 1:n
			for j = 1:i
                innerSum = getInnerSum(i, VGH(1), x0, innerSums);
				outGr(j) = outGr(j) + 2*innerSum;
                innerSums(i) = innerSum;
			end
		end
    end
    
    if VGH(3) > 0
		outHes = zeros(n,n);

        for i = n:-1:1
			for j = i:-1:1
                for k = j:n
                    outHes(n-i+1,n-k+1) = outHes(n-i+1,n-k+1) + 2;
                end
			end
        end
    else
        outHes = 0;
    end;
    
end
    

function is = getInnerSum(i, alreadyComputed, x, innerSums)
    if alreadyComputed
        is = innerSums(i);
    else
        is = -i;
        is = is + sum(x(1:i));
    end
end

