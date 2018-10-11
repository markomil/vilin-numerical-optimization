function [ outVal,outGr,outHes ] = FullHessian1( x0,VGH )
% Full Hessian FH1 function

	n = length(x0);
	outVal = 0;

	% computes the value of function in point x0
	if VGH(1) > 0
		outVal = (x0(1)-3).^2;
		sum = x0(1);
		for i=2:n
			sum = sum + x0(i);
			outVal = outVal + (x0(1)-3-2*sum.^2).^2;
		end
	end

	% computes the numerical gradient value of function in point x0
	if VGH(2) > 0
	   outGr = zeros(n, 1);
	   pref = zeros(n,1);
	   postWithCoefs = zeros(n,1);
	   pref(1) = x0(1);
	   postWithCoefs(n) = x0(n);
	   for i=2:n
		  pref(i) = pref(i-1) + x0(i);
	   end
	   postCubeOfPref = zeros(n,1);
	   postCubeOfPref(n) = pref(n).^3;
	   for i=n-1:-1:1
		  postCubeOfPref(i) = postCubeOfPref(i+1) + pref(i).^3;
		  postWithCoefs(i) = postWithCoefs(i+1)+(n-i+1)*x0(i); 
	   end
	   for i=2:n-1
		  outGr(i) = -8*((x0(1)-3)*((n-i+1)*pref(i)+postWithCoefs(i+1))-2*postCubeOfPref(i)); 
	   end
	   outGr(1) = 2*(x0(1)-3);
	   for i=2:n
		  outGr(1) = outGr(1) + 2*(x0(1)-3-2*pref(i).^2)*(1-4*pref(i));
	   end
	   outGr(n) = -8*((x0(1)-3)*pref(n)-2*postCubeOfPref(n));
	else
	   outGr = 0;
	end

	% computes the numerical Hessian of function in point x0
	if VGH(3) > 0
	   outHes = zeros(n, n);     
	   postQuadPref = zeros(n,1);
	   pref = zeros(n,1);
	   pref(1) = x0(1);
	   for i=2:n
		  pref(i) = pref(i-1) + x0(i);
	   end
	   postQuadPref(n) = pref(n).^2;
	   for i=n-1:-1:1
		   postQuadPref(i) = postQuadPref(i+1)+pref(i).^2;
	   end
	   for i=2:n-1
		   for j=i:n-1
			   outHes(i,j) = -8*(n-j+1)*(x0(1)-3)+48*postQuadPref(j);
			   outHes(j,i) = outHes(i,j);
		   end
	   end
	   
	   postSumPref=zeros(n,1);
	   postSumPref(n) = pref(n);
	   for i=n-1:-1:1
		  postSumPref(i)= postSumPref(i+1)+pref(i); 
	   end
	   
	   for i=2:n
		   outHes(1,i) = outHes(i,i)-8*postSumPref(i);
		   outHes(i,1) = outHes(1,i);
	   end
	   for i=2:n
		   outHes(n,i) = -8*(x0(1)-3 - 6*pref(n).^2);
		   outHes(i,n) = outHes(n,i);
	   end
	   outHes(1,n) = -8*(x0(1)-3+pref(n)-6*pref(n).^2);
	   outHes(n,1) = outHes(1,n);
	   
	   outHes(1,1)=2;
	   for i=2:n
		   outHes(1,1) = outHes(1,1)+2*(13-8*pref(i)+24*pref(i).^2-4*x0(1));
	   end
	else
		outHes = 0;
	end

end

