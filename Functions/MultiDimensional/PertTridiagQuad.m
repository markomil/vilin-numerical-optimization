function [ outVal, outGr, outHes ] = PertTridiagQuad( x0, VGH )
% Perturbed Tridiagonal Quadratic function

	n = length(x0);
	outVal = 0;
	outGr = 0;
	
	if VGH(1) > 0
		outVal = x0(1)^2;
		for i=2:n-1
			outVal = outVal + i*x0(i)^2 + (x0(i-1)+x0(i)+x0(i+1))^2;
		end
	end
	
	if VGH(2) > 0
		outGr = zeros(n,1);
		outGr(1) = 2*x0(1);
		for i = 2:n-1
			outGr(i-1) = outGr(i-1) + 2*(x0(i-1)+x0(i)+x0(i+1));
			outGr(i) = outGr(i) + 2*i*x0(i) + 2*(x0(i-1)+x0(i)+x0(i+1));
			outGr(i+1) = outGr(i+1) + 2*(x0(i-1)+x0(i)+x0(i+1));
		end
    end

    if VGH(3) > 0
		outHes = zeros(n,n);

        % first row
		outHes(1,1) = 4;
        outHes(1,2) = 2;
        outHes(1,3) = 2;
        
        % second row
		outHes(2,1) = 2;
        outHes(2,2) = 8;
        outHes(2,3) = 4;
        outHes(2,4) = 2;
                
        for i = 3:n-2
			outHes(i,i-2) = 2;
			outHes(i,i-1) = 4;
			outHes(i,i) = 2*i+6;
            
			outHes(i,i+1) = 4;
			outHes(i,i+2) = 2;
        end
        
        % row n-1
        outHes(n-1,n-3) = 2;
        outHes(n-1,n-2) = 4;
        outHes(n-1,n-1) = 2*(n-2)+6;
        outHes(n-1,n) = 2;
        
        % n-th row
        outHes(n,n-2) = 2;
        outHes(n,n-1) = 2;
        outHes(n,n) = 2;
    else
        outHes = 0;
    end;
    
end
