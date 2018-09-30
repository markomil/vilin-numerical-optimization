function [ outVal,outGr,outHes ] = mojaExtFreudAndRoth( x0,VGH )

n = length(x0);

assert (mod(n,2)==0)

outVal = 0;
outGr = zeros(n, 1);

if VGH(1) > 0
   for i=2:2:n
       outVal = outVal + (-13 + x0(i-1) + ((5-x0(i))*x0(i)-2)*x0(i)).^2 + ...
           (-29+x0(i-1)+((x0(i)+1)*x0(i)-14)*x0(i)).^2;
   end
end

if VGH(2) > 0
   for i=1:n
      if(mod(i,2) == 1)
          outGr(i) = 4*x0(i)+12*x0(i+1).^2-32*x0(i+1)-84;
      else
          outGr(i) = 24*x0(i-1)*x0(i)-32*x0(i-1)+12*x0(i).^5- ...
              40*x0(i).^4+8*x0(i).^3-240*x0(i).^2+24*x0(i)+864;
      end
   end
end

if VGH(3) > 0
   outHes = zeros(n, n);     
   for i=2:2:n
       outHes(i-1,i-1) = 4;
       outHes(i,i) = 24*x0(i-1)+60*x0(i).^4-160*x0(i).^3 + ...
                    24*x0(i).^2-480*x0(i)+24;
       outHes(i-1,i) = 24*x0(i)-32;
       outHes(i,i-1) = 24*x0(i)-32;
   end
else
    outHes = 0;
end

end

