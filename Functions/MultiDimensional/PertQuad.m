function [ outVal, outGr, outHes ] = PertQuad( x, VGH )
% Perturbed Quadratic function

    c = 100;
    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
    sumX = sum(x);
    
    if VGH(1) > 0
        for i = 1:n
            outVal = outVal + i*(x(i)^2);
        end
        outVal = outVal + 1/c*sumX^2;
    end

    if VGH(2) > 0
        for i = 1:n
            outGr(i) = 2*i*x(i) + 2/c*sumX;
        end
    end

    if VGH(3) > 0
        outHes = 2/c*ones(n, n);
        for i = 1:n
           outHes(i, i) = 2*i + 2/c;
        end
    else
        outHes = 0;
    end;

end
