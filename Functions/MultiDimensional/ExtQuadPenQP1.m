function [ outVal, outGr, outHes ] = ExtQuadPenQP1( x, VGH )
% Extended quadratic penalty QP1 function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    exp2 = sum(x.^2)-0.5;
    
    if VGH(1) > 0
        exp1 = sum((x(1:n-1).^2 - 2).^2);
        outVal = exp1 + exp2^2;
    end

    if VGH(2) > 0
        for i = 1:n-1
            outGr(i) = 4*x(i)*(-2 + x(i)^2) + 4*x(i)*exp2;
        end
        outGr(n) = 4*x(n)*exp2;
    end
          
    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            for j = 1:n
                if i == j
                    outHes(i, j) =  16*x(i)^2 + 4*(-2 + x(i)^2) + 4*exp2;
                else
                    outHes(i, j) = 8*x(i)*x(j);
                end
            end
        end;
        outHes(n, n) = 8*x(n)^2 + 4*exp2;
    else
        outHes = 0;
    end

end
