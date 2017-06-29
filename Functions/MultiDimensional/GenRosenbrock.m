function [outVal, outGr, outHes ] = GenRosenbrock( x0, VGH)
% Generalized Rosenbrock function

    dim = length(x0);    
    coef = 100;
    outVal = 0;
    outGr = zeros(dim,1);


    if VGH(1) >0
        for i=1:dim-1
            outVal = outVal + (1-x0(i)).^2 + coef*(x0(i+1) - x0(i).^2).^2;
        end
    end

    if VGH(2) >0
        outGr(1) = outGr(1) - 2*(1-x0(1)) + coef*(4*x0(1).^3 - 4*x0(1)*x0(2));
        for i=2:dim-1
            outGr(i) = outGr(i) - 2*(1-x0(i)) + coef*(4*x0(i).^3 - 2*x0(i-1).^2 - 4*x0(i)*x0(i+1) + 2*x0(i));
        end
        outGr(dim) = outGr(dim) + coef*(2*x0(dim) - 2*x0(dim-1).^2);
    end

    % computes the numerical values of Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(dim,dim);
        % computes first row
        outHes(1,1) =  double(2 + 8*coef*x0(1)^2 - 4*coef*(-x0(1)^2 + x0(2)));
        outHes(1,2) = -4*coef*x0(1);
        % computes rows from 2 to dim-1
        for i=2:dim-1
            outHes(i,i-1) = double(-4*coef*x0(i-1));
            outHes(i,i) = double(2 + 2*coef + 8*coef*x0(i)^2 - 4*coef*(-x0(i)^2 + x0(i+1)));
            outHes(i,i+1) = double(-4*coef*x0(i));
        end
        % computes last row 
        outHes(dim,dim) = double(2*coef);
        outHes(dim,dim-1) = -4*coef*x0(dim-1);
    else
        outHes = 0;
    end;
    
end
