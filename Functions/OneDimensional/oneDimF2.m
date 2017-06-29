function [ outVal, outGr ] = oneDimF2( x )
%Function f2 from functions.txt

outVal = -x^3*exp(-x);
outGr = x^2*exp(-x)*(x-3);

end

