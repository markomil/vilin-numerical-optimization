function [ outVal, outGr, outHes ] = SINE( x, VGH )

    n = length(x);
    outVal = 0;
    outGr = zeros(n, 1);
	
    if VGH(1) > 0
        for i = 1:n-1
            outVal = outVal + sin(-0.5*x(i+1) + x(i)^2);
        end
        
    end

    if VGH(2) > 0
        for i = 2:n-1
            outGr(i) = 2*x(i)*cos(-0.5*x(i+1) + x(i)^2) - 0.5*cos(-0.5*x(i) + x(i-1)^2);
        end
        outGr(1) = 2*x(1)*cos(-0.5*x(2) + x(1)^2);
        outGr(n) = -0.5*cos(-0.5*x(n) + x(n-1)^2);
    end

    if VGH(3) > 0
        outHes = zeros(n, n);
        outHes(1,1) =  2*cos(-0.5*x(2) + x(1)^2)-4*x(1)^2*sin(-0.5*x(2) + x(1)^2);
        outHes(1,2) = 2*0.5*sin(-0.5*x(2) + x(1)^2);
        outHes(2,1) = outHes(1,2);
        outHes(n,n) = -0.5*0.5*sin(-0.5*x(n) + x(n-1)^2);
        %outHes(n,n-1) = -2*0.5*x(n-1)*sin(-0.5*x(n) + x(n-1)^2);
        for i=2:n-1
            outHes(i+1,i) = 2*x(i)*0.5*sin(-0.5*x(i+1) + x(i)^2);
            outHes(i,i) = 2*cos(-0.5*x(i+1) + x(i)^2) -4*x(i)^2*sin(-0.5*x(i+1) + x(i)^2) -0.25*sin(-0.5*x(i) + x(i-1)^2);
            outHes(i,i+1) = outHes(i+1,i);
        end
    else
        outHes = 0; 
    end
end
