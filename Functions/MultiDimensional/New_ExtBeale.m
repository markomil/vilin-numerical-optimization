function [ outVal,outGr,outHes ] = mojaExtBeale( x0,VGH )

n = length(x0);

assert (mod(n,2)==0)

outVal = 0;
outGr = zeros(n, 1);

if VGH(1) > 0
   for i=2:2:n
       outVal = outVal + (1.5-x0(i-1)*(1-x0(i))).^2 + (2.25-x0(i-1)*(1-x0(i).^2)).^2 + ...
           (2.625-x0(i-1)*(1-x0(i).^3)).^2;
   end
end

if VGH(2) > 0
   for i=1:n
      if(mod(i,2) == 1)
          outGr(i) = 2*(x0(i+1)-1)*(x0(i)*(x0(i+1)-1)+1.5) + 2*(x0(i+1).^2 - 1)*(x0(i)*(x0(i+1).^2 - 1) + 2.25) + ...
              2*(x0(i+1).^3 - 1)*(x0(i)*(x0(i+1).^3 - 1) + 2.625);
      else
          outGr(i) = 2*x0(i-1)*(x0(i-1)*(x0(i) - 1) + 3*x0(i).^2*(x0(i-1)*(x0(i).^3 - 1) + 2.625) + 2*x0(i)*(x0(i-1)*(x0(i).^2 - 1) + 2.25) + 1.5);
      end
   end
end

if VGH(3) > 0
   outHes = zeros(n, n);     
   for i=2:2:n
       outHes(i-1,i-1) = 2*(x0(i) - 1).^2 + 2*(x0(i).^2 - 1).^2 + 2*(x0(i).^3 - 1).^2;
       outHes(i,i-1) = 12.0*x0(i-1)*x0(i).^5 + 8.0*x0(i-1)*x0(i).^3 - 12.0*x0(i-1)*x0(i).^2 - 4.0*x0(i-1)*x0(i) - 4.0*x0(i-1) + 15.75*x0(i).^2 + 9.0*x0(i) + 3.0;
       outHes(i-1,i) = 12.0*x0(i-1)*x0(i).^5 + 8.0*x0(i-1)*x0(i).^3 - 12.0*x0(i-1)*x0(i).^2 - 4.0*x0(i-1)*x0(i) - 4.0*x0(i-1) + 15.75*x0(i).^2 + 9.0*x0(i) + 3.0;
       outHes(i,i) = x0(i-1)*(30.0*x0(i-1)*x0(i).^4 + 12.0*x0(i-1)*x0(i).^2 - 12.0*x0(i-1)*x0(i) - 2.0*x0(i-1) + 31.5*x0(i) + 9.0);
   end
else
    outHes = 0;
end


end

