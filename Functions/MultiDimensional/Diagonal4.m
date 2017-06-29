function [ outVal, outGr, outHes ] = Diagonal4( x, VGH )
% Diagonal 4 function

    c = 100;
    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    
    if VGH(1) > 0
        for i = 1:n/2
            outVal = outVal + x(2*i-1)^2 + c*x(2*i)^2;
        end
        outVal = outVal/2;
    end

    if VGH(2) > 0
        for i = 1:n/2
            outGr(2*i-1) = x(2*i-1);
            outGr(2*i) = c * x(2*i);
        end
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        for i = 1:n/2
            outHes(2*i-1, 2*i-1) = 1;
            outHes(2*i, 2*i) = c;
        end
    else
        outHes = 0;
    end

end
