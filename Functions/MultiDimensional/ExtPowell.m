function [ outVal, outGr, outHes ] = ExtPowell( x0, VGH )
% Extended Powell function

    n = length(x0);
    outVal = 0;
    outGr = 0;
    
    if VGH(1) > 0
    	for i=1:n/4
			outVal = outVal + (x0(4*i-3)+10*x0(4*i-2))^2 + 5*(x0(4*i-1)-x0(4*i))^2 + (x0(4*i-2)-2*x0(4*i-1))^4 + 10*(x0(4*i-3)-x0(4*i))^4;
    	end
    end


    if VGH(2) > 0
    		outGr = zeros(n, 1);
		for i=1:n/4
			outGr(4*i-3) = 2*(x0(4*i-3)+10*x0(4*i-2)) + 40*(x0(4*i-3)-x0(4*i))^3;
			outGr(4*i-2) = 20*(x0(4*i-3)+10*x0(4*i-2)) + 4*(x0(4*i-2)-2*x0(4*i-1))^3;
			outGr(4*i-1) = 10*(x0(4*i-1)-x0(4*i)) - 8*(x0(4*i-2)-2*x0(4*i-1))^3;
			outGr(4*i) = -10*(x0(4*i-1)-x0(4*i)) - 40*(x0(4*i-3)-x0(4*i))^3;
		end	
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
		for i=1:n/4
			outHes(4*i-3, 4*i-3) = 2 + 120*(x0(4*i-3)-x0(4*i))^2;
			outHes(4*i-3, 4*i-2) = 20;
			outHes(4*i-3, 4*i) = -120*(x0(4*i-3)-x0(4*i))^2;

			outHes(4*i-2, 4*i-3) = 20;
			outHes(4*i-2, 4*i-2) = 200 + 12*(x0(4*i-2)-2*x0(4*i-1))^2;
			outHes(4*i-2, 4*i-1) = -24*(x0(4*i-2)-2*x0(4*i-1))^2;

			outHes(4*i-1, 4*i-2) = -24*(x0(4*i-2)-2*x0(4*i-1))^2;
			outHes(4*i-1, 4*i-1) = 10 + 48*(x0(4*i-2)-2*x0(4*i-1))^2;
			outHes(4*i-1, 4*i) = -10;

			outHes(4*i, 4*i-3) = -120*(x0(4*i-3)-x0(4*i))^2;
			outHes(4*i, 4*i-1) = -10;
			outHes(4*i, 4*i) = 10 + 120*(x0(4*i-3)-x0(4*i))^2;
		end
    else
        outHes = 0;
    end;

end

