function [ outVal,outGr,outHes ] = mojaExtTET( x0,VGH )

n = length(x0);

assert (mod(n,2)==0)

outVal = 0;
outGr = zeros(n, 1);

if VGH(1) > 0
   for i=2:2:n
       outVal = outVal + exp(x0(i-1)+3*x0(i)-0.1)+exp(x0(i-1)-3*x0(i)-0.1)+exp(-x0(i-1)-0.1);
   end
end

if VGH(2) > 0
   for i=1:n
      if mod(i,2) == 1
          outGr(i) = 0.90483741803596*((exp(x0(i) - 3*x0(i+1)) + exp(x0(i) + 3*x0(i+1)))*exp(x0(i)) - 1)*exp(-x0(i));
      else
          outGr(i) = -2.71451225410788*exp(x0(i-1) - 3*x0(i)) + 2.71451225410788*exp(x0(i-1) + 3*x0(i));
      end
   end
end

if VGH(3) > 0
   outHes = zeros(n, n);     
   for i=2:2:n
       outHes(i-1,i-1) = 0.90483741803596*((exp(x0(i-1) - 3*x0(i)) + exp(x0(i-1) + 3*x0(i)))*exp(x0(i-1)) + 1)*exp(-x0(i-1));
       outHes(i,i)=8.14353676232364*exp(x0(i-1) - 3*x0(i)) + 8.14353676232364*exp(x0(i-1) + 3*x0(i));
       outHes(i,i-1) = -2.71451225410788*exp(x0(i-1) - 3*x0(i)) + 2.71451225410788*exp(x0(i-1) + 3*x0(i));
       outHes(i-1,i) = -2.71451225410788*exp(x0(i-1) - 3*x0(i)) + 2.71451225410788*exp(x0(i-1) + 3*x0(i));
   end
else
    outHes = 0;
end


end

