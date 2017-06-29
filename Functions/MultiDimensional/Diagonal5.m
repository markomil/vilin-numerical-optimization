function [ outVal, outGr, outHes ] = Diagonal5( x, VGH )
% Diagonal 5 function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + log((exp(x(i)) + exp(-x(i))));
        end
    end

    if VGH(2) > 0
        for i = 1:n
            ex1 = exp(2*x(i));
            outGr(i) = (ex1 - 1)/(ex1 + 1);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            ex1 = exp(2*x(i));
            outHes(i, i) = 4*ex1/(ex1 + 1)^2;
        end
    else
        outHes = 0;
    end

end

