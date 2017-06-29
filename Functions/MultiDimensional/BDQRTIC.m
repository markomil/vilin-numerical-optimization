function [ outVal, outGr, outHes ] = BDQRTIC( x0, VGH )
% BDQRTIC (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
    if VGH(1) > 0
        for i=1:n-4
			outVal = outVal + (-4*x0(i) + 3)^2 + (x0(i)^2 + 2*x0(i+1)^2 + 3*x0(i+2)^2 + 4*x0(i+3)^2 + 5*x0(n)^2)^2;
        end
    end
    
    if VGH(2) > 0
		outGr = zeros(n,1);

		ba1 = bigAddend(1, x0, n);
		ba2 = bigAddend(2, x0, n);
		ba3 = bigAddend(3, x0, n);
		ba4 = bigAddend(4, x0, n);

		outGr(1) = -8*(-4*x0(1) + 3) + 4*x0(1)*ba1;
		outGr(2) = 8*x0(2)*ba1 - 8*(-4*x0(2) + 3) + 4*x0(2)*ba2;
		outGr(3) = 12*x0(3)*ba1 + 8*x0(3)*ba2 + 4*x0(3)*ba3 - 8*(-4*x0(3) + 3);
		outGr(n) = 20*x0(n)*ba1 + 20*x0(n)*ba2 + 20*x0(n)*ba3;

        for i = 4:n-4
			outGr(i) = 16*x0(i)*ba1 + 12*x0(i)*ba2 + 8*x0(i)*ba3 + 4*x0(i)*ba4 - 8*(-4*x0(i) + 3);
			outGr(n) = outGr(n) + 20*x0(n)*ba4;

			ba1 = ba2;
			ba2 = ba3;
			ba3 = ba4;
			ba4 = bigAddend(i+1, x0, n);
        end
        
        outGr(n-3) = 16*x0(n-3)*ba1 + 12*x0(n-3)*ba2 + 8*x0(n-3)*ba3;
		outGr(n-2) = 16*x0(n-2)*ba2 + 12*x0(n-2)*ba3;
		outGr(n-1) = 16*x0(n-1)*ba3;
    end
    
    if VGH(3) > 0
		outHes = zeros(n,n);

		ba1 = bigAddend(1, x0, n);
		ba2 = bigAddend(2, x0, n);
		ba3 = bigAddend(3, x0, n);
		ba4 = bigAddend(4, x0, n);

		outHes(1,1) = 32 + 8*x0(1)^2 + 4*ba1;
		outHes(1,2) = 16*x0(1)*x0(2);
		outHes(1,3) = 24*x0(1)*x0(3);
		outHes(1,4) = 32*x0(1)*x0(4);
		outHes(1,n) = outHes(1,n) + 40*x0(1)*x0(n);

		outHes(2,1) = 16*x0(1)*x0(2);
		outHes(2,2) = 32 + 40*x0(2)^2 + 8*ba1 + 4*ba2;
		outHes(2,3) = 64*x0(2)*x0(3);
		outHes(2,4) = 88*x0(2)*x0(4);
		outHes(2,5) = 32*x0(2)*x0(5);
		outHes(2,n) = outHes(2,n) + 120*x0(2)*x0(n);

		outHes(3,1) = 24*x0(3)*x0(1);
		outHes(3,2) = 64*x0(2)*x0(3);
		outHes(3,3) = 112*x0(3)^2 + 12*ba1 + 8*ba2 + 4*ba3 + 32;
		outHes(3,4) = 160*x0(3)*x0(4);
		outHes(3,5) = 88*x0(3)*x0(5);
		outHes(3,6) = 32*x0(3)*x0(6);
		outHes(3,n) = outHes(3,n) + 240*x0(3)*x0(n);

		outHes(n,1) = 40*x0(1)*x0(n);
		outHes(n,2) = 120*x0(2)*x0(n);
		outHes(n,3) = 240*x0(3)*x0(n);
		outHes(n,4) = 360*x0(4)*x0(n);
		outHes(n,5) = 280*x0(5)*x0(n);
		outHes(n,6) = 160*x0(6)*x0(n);
		outHes(n,n) = 20*(ba1 + ba2 + ba3) + 600*x0(n)^2;

		for i = 4:n-4
			outHes(i,i-3) = 32*x0(i)*x0(i-3);
			outHes(i,i-2) = 88*x0(i)*x0(i-2);
			outHes(i,i-1) = 160*x0(i)*x0(i-1);
			outHes(i,i)   = 240*x0(i)^2 + 16*ba1 + 12*ba2 + 8*ba3 + 4*ba4 + 32;
			outHes(i,i+1) = 160*x0(i)*x0(i+1);
			outHes(i,i+2) = 88*x0(i)*x0(i+2);
			outHes(i,i+3) = 32*x0(i)*x0(i+3);
			outHes(i,n)   = outHes(i,n) + 400*x0(i)*x0(n);

			outHes(n,i)   = outHes(n,i) + 40*x0(i)*x0(n);
			outHes(n,i+1) = outHes(n,i+1) + 80*x0(i+1)*x0(n);	
			outHes(n,i+2) = outHes(n,i+2) + 120*x0(i+2)*x0(n);
			outHes(n,i+3) = outHes(n,i+3) + 160*x0(i+3)*x0(n);
			outHes(n,n)   = outHes(n,n) + 20*ba4 + 200*x0(n)^2;

			ba1 = ba2;
			ba2 = ba3;
			ba3 = ba4;
			ba4 = bigAddend(i+1, x0, n);
		end

		outHes(n-3, n-6) = outHes(n-3, n-6) + 32*x0(n-3)*x0(n-6);
		outHes(n-3, n-5) = outHes(n-3, n-5) + 88*x0(n-3)*x0(n-5);
		outHes(n-3, n-4) = outHes(n-3, n-4) + 160*x0(n-3)*x0(n-4);
		%outHes(n-3, n-3) = outHes(n-3, n-3) + 232*x0(n-3)^2 + 16*ba2 + 12*ba3 + 8*ba4;
        outHes(n-3, n-3) = outHes(n-3, n-3) + 232*x0(n-3)^2 + 16*ba1 + 12*ba2 + 8*ba3;
		outHes(n-3, n-2) = outHes(n-3, n-2) + 144*x0(n-3)*x0(n-2);
		outHes(n-3, n-1) = outHes(n-3, n-1) + 64*x0(n-3)*x0(n-1);
		outHes(n-3, n)   = outHes(n-3, n) + 360*x0(n-3)*x0(n);

		outHes(n-2, n-5) = outHes(n-2, n-5) + 32*x0(n-2)*x0(n-5);
		outHes(n-2, n-4) = outHes(n-2, n-4) + 88*x0(n-2)*x0(n-4);
		outHes(n-2, n-3) = outHes(n-2, n-3) + 144*x0(n-2)*x0(n-3);
		%outHes(n-2, n-2) = outHes(n-2, n-2) + 200*x0(n-2)^2 + 16*ba3 + 12*ba4;
        outHes(n-2, n-2) = outHes(n-2, n-2) + 200*x0(n-2)^2 + 16*ba2 + 12*ba3;
		outHes(n-2, n-1) = outHes(n-2, n-1) + 96*x0(n-2)*x0(n-1);
		outHes(n-2, n)   = outHes(n-2, n) + 280*x0(n-2)*x0(n);

		outHes(n-1, n-4) = outHes(n-1, n-4) + 32*x0(n-1)*x0(n-4);
		outHes(n-1, n-3) = outHes(n-1, n-3) + 64*x0(n-1)*x0(n-3);
		outHes(n-1, n-2) = outHes(n-1, n-2) + 96*x0(n-1)*x0(n-2);
		%outHes(n-1, n-1) = outHes(n-1, n-1) + 128*x0(n-1)^2 + 16*ba4;
        outHes(n-1, n-1) = outHes(n-1, n-1) + 128*x0(n-1)^2 + 16*ba3;
		outHes(n-1, n)   = outHes(n-1, n) + 160*x0(n-1)*x0(n);
    else
        outHes = 0;
    end;
    
end


function [ val ] = bigAddend(i, x0, n)
% Big addend of the function
    val = (x0(i)^2 + 2*x0(i+1)^2 + 3*x0(i+2)^2 + 4*x0(i+3)^2 + 5*x0(n)^2);
end
