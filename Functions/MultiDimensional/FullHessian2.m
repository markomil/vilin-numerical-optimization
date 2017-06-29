function [ outVal, outGr, outHes ] = FullHessian2( x0, VGH )
% Full Hessian 2 function
    
	n=length(x0);
	outVal = 0;
	outGr = 0;
	
	% constants
	c1 = 5;
	c2 = 1;
       
	if VGH(1) > 0
		outVal = (x0(1)-c1)^2;
		innerVal = x0(1)-c2;
		for i=2:n
			innerVal = innerVal + x0(i);
			outVal = outVal + innerVal^2;
		end	
	end
    
	if VGH(2) > 0
		outGr = zeros(n, 1);

		oldInnerVal = 0;
		innerVal = sum(x0) - c2;
		for i=n:-1:2
			outGr(i) = 2*(innerVal+oldInnerVal);
			oldInnerVal = oldInnerVal + innerVal;
			innerVal = innerVal - x0(i);
		end
		outGr(1) = 2*(x0(1) - c1) + 2*oldInnerVal;
    end
    
    if VGH(3) > 0
		outHes = zeros(n, n);
		row = repmat(2, 1, n);
		tmpRow = row;
        for i=n:-1:1
			outHes(i, :) = row;
			tmpRow(i) = 0;
			row = row + tmpRow;
        end
    else
        outHes = 0;
    end;

end

