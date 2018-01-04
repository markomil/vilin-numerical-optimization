function [outVal, outGr, outHes ] = SINQUADNew(x0, VGH)
% SINQUAD (CUTE) function
% TODO check if optimal computation of gradient and Hessian is done

    dim = length(x0);    
    outVal = 0;
    outGr = 0;
    outHes = 0;

    % computes the value of function in point x0
    if VGH(1) > 0
        outVal = (x0(1)-1)^4+(x0(dim)^2-x0(1)^2)^2;
        
        for i = 2:dim-1
            outVal = outVal + (sin(x0(i)-x0(dim))-x0(1)^2+x0(i)^2)^2;
        end
    end
    
    % computes the numerical gradient value of function in point x0
    if VGH(2) > 0
        outGr = zeros(dim,1);
        outGr(1) = 4*((-1 + x0(1))^3 + x0(1)^3 - x0(1)*x0(dim)^2);
        outGr(dim)=4*x0(dim)*(x0(dim)^2-x0(1)^2);
        
        for i = 2:dim-1
            outGr(1) = outGr(1)-4*x0(1)*(-x0(1)^2 + x0(i)^2 + sin(x0(i) - x0(dim)));
            outGr(i) = 2*(cos(x0(i) - x0(dim)) + 2*x0(i))*(sin(x0(i)-x0(dim)) - x0(1)^2 + x0(i)^2);
            outGr(dim) = outGr(dim) -2*cos(x0(i)-x0(dim))*(-x0(1)^2+x0(i)^2+sin(x0(i)-x0(dim)));
        end
    
    end

    % computes the numerical Hessian of function in point x0
    if VGH(3) > 0
        outHes = zeros(dim,dim);
        
        % racunanje vrendosti Hesijana na poziciji (1,1) 
        outHes(1,1) =  4*(3 + 6*(-1 + x0(1))*x0(1) - x0(dim)^2);
        
        for i = 2:dim-1
            outHes(1,1) = outHes(1,1) +12*x0(1)^2 - 4*x0(i)^2 - 4*sin(x0(i) - x0(dim));
        end

        % racunanje vrendosti Hesijana na poziciji (dim,dim) 
        outHes(dim,dim) =  -4*(x0(1)^2 - 3*x0(dim)^2);
        
        for i = 2:dim-1
            outHes(dim,dim) = outHes(dim,dim) + 2*(cos(2*(x0(i) - x0(dim))) + (x0(1) - x0(i))*(x0(1) + x0(i))*sin(x0(i) - x0(dim)));
        end

        outHes(1,dim) =  -8*x0(1)*x0(dim);
        
        for i = 2:dim-1
            outHes(1,dim) = outHes(1,dim) + 4*x0(1)*cos(x0(i)-x0(dim)) ;
        end
        
        outHes(dim,1)=outHes(1,dim);

        % racunanje vrendosti Hesijana na pozicijama (1,i)
        for i = 2:dim-1        
            outHes(1,i) =-4*x0(1)*(2*x0(i)+cos(x0(i)-x0(dim)));
            outHes(i,1) = outHes(1,i);        
        end

        % racunanje vrendosti Hesijana na pozicijama (dim,i) i (i,dim)
        for i = 2:dim-1
            outHes(dim,i) =-2*cos(x0(i)-x0(dim))*(2*x0(i)+cos(x0(i)-x0(dim)))+2*sin(x0(i)-x0(dim))*(-x0(1)^2+x0(i)^2+sin(x0(i)-x0(dim)));
            outHes(i,dim) = outHes(dim,i);        
        end

        % racunanje vrendosti Hesijana na pozicijama (i,j) za i,j=2,...,dim-1
        for i = 2:dim-1         
            outHes(i,i) =2*(2*x0(i)+cos(x0(i)-x0(dim)))^2+2*(2-sin(x0(i)-x0(dim)))*(-x0(1)^2+x0(i)^2+sin(x0(i)-x0(dim)));      
        end
    end
    
end

    


