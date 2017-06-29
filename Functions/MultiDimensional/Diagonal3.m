function [ outVal, outGr, outHes ] = Diagonal3( x, VGH )
% Diagonal 3 function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + (exp(x(i)) - i*sin(x(i)));
        end
    end

    if VGH(2) > 0
        for i = 1:n
            outGr(i) = exp(x(i)) - i*cos(x(i));
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            outHes(i, i) = exp(x(i)) + i*sin(x(i));
        end
    else outHes = 0;
    end

end
