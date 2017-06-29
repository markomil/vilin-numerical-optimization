function [ outVal, outGr, outHes ] = ExtPen( x0, VGH )
% Extended Penalty function
% Check the initial point x0 = (1, 2, ... n)?!

    n = length(x0);
    outVal = 0;
    outGr = zeros(n, 1);
    sum1 = sum((x0(1:n-1)-1).^2);
    sum2 = sum(x0.^2-0.25);
    
    if VGH(1) > 0
        outVal = sum1 + sum2^2;
    end

    if VGH(2) > 0
        for i = 1:n-1
            outGr(i) = 2*(x0(i)-1) + 4*sum2*x0(i);
        end
        outGr(n) = 4*sum2*x0(n);
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n-1
           outHes(i, i) = 2 + 4*sum2 + 8*x0(i)^2;
           for j = i+1:n
               outHes(i, j) = 8*x0(i)*x0(j);
               outHes(j, i) = 8*x0(i)*x0(j);
           end
           outHes(n, n) = 4*sum2 + 8*x0(n)^2;
        end
    else
        outHes = 0;
    end;

end
