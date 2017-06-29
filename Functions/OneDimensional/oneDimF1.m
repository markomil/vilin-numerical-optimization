function [ outVal, outGr ] = oneDimF1( x )
%Function f1 from functions.txt

outVal = -x^2*exp(-x);
outGr = x*exp(-x)*(x-2);

end

