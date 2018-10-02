function [ outVal,outGr,outHes ] = New_Staircase1( x0,VGH )
% Staircase 1 function

	n = length(x0);
	outVal = 0;
	outGr = zeros(n, 1);
	sum = 0;

	% computes the value of function in point x0
	if VGH(1) > 0
		for i=1:n
			sum = sum + x0(i);
			outVal = outVal + sum.^2;
		end
	end

	% computes the numerical gradient value of function in point x0
	if VGH(2) > 0
	   pref = zeros(n,1);
	   postWithCoefs = zeros(n,1);
	   pref(1) = x0(1);
	   postWithCoefs(n) = x0(n);
	   for i=2:n
		  pref(i) = pref(i-1) + x0(i);
	   end
	   for i=n-1:-1:1
		  postWithCoefs(i) = postWithCoefs(i+1)+(n-i+1)*x0(i); 
	   end
	   for i=2:n-1
		  outGr(i) = 2*((n-i+1)*pref(i)+postWithCoefs(i+1)); 
	   end
	   outGr(1) = 2*postWithCoefs(1);
	   outGr(n) = 2*pref(n);
	end

	% computes the numerical Hessian of function in point x0
	if VGH(3) > 0
	   outHes = zeros(n, n);
	   for i=1:n
		   outHes(i,i) = 2*(n-i+1);
		   for j=i+1:n
			   outHes(i,j) = 2*(n-j+1);
			   outHes(j,i) = 2*(n-j+1);
		   end
	   end
	else
		outHes = 0;
	end

end

