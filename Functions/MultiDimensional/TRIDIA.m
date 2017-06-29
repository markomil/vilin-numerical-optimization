function [ outVal, outGr, outHes ] = TRIDIA( x0, VGH )
% TRIDIA (CUTE) function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	% constants
	alpha = 2;
	beta = 1;
	gama = 1;
	delta = 1;

	if VGH(1) > 0
		outVal = gama*(delta*x0(1)-1)^2;
		for i = 2:n
			outVal = outVal + i*(alpha*x0(i) - beta*x0(i-1))^2;
		end
    end
	
    if VGH(2) > 0
		outGr = zeros(n,1);
		outGr(1) = 2*delta*gama*(delta*x0(1)-1) - 2*2*beta*(alpha*x0(2) - beta*x0(1));
        
        for i = 2:n-1
			outGr(i) = 2*i*alpha*(alpha*x0(i) - beta*x0(i-1)) - 2*(i+1)*beta*(alpha*x0(i+1) - beta*x0(i));
        end
        outGr(n) = 2*n*alpha*(alpha*x0(n) - beta*x0(n-1));
    end

    if VGH(3) > 0
		outHes = zeros(n,n);
		outHes(1,1) = 2*gama*delta^2 + 4*beta^2;
		outHes(1,2) = -4*alpha*beta;
        
        for i = 2:n-1
			outHes(i,i-1) = -2*i*alpha*beta;
			outHes(i,i)   = 2*i*alpha^2 + 2*(i+1)*beta^2;
			outHes(i,i+1) = -2*(i+1)*alpha*beta;
        end
        
		outHes(n,n-1) = -2*n*alpha*beta;
		outHes(n,n)   = 2*n*alpha^2;
    else
        outHes = 0;
    end;
       
end

