function [ outVal, outGr, outHes ] = Hager( x, VGH )
% Hager function

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + exp(x(i)) - sqrt(i)*x(i);
        end
    end

    if VGH(2) > 0
        for i = 1:n
            outGr(i) = exp(x(i)) - sqrt(i);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            outHes(i,i) = exp(x(i));
        end
    else
        outHes = 0;
    end

end
