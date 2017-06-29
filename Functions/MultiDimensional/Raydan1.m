function [ outVal, outGr, outHes ] = Raydan1( x, VGH )
% Raydan 1 function

    param = 10;
    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);

    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + (i/param)*(exp(x(i)) - x(i));
        end
    end

    if VGH(2) > 0
        for i = 1:n
            outGr(i) = (i/param)*(exp(x(i)) - 1);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n
            outHes(i,i) = i/10*(exp(x(i)));
        end
    else
        outHes = 0;
    end;

end
