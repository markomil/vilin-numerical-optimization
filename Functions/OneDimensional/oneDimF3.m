function [ outVal, outGr ] = oneDimF3( x )
%Function f3 from functions.txt

outVal = -x^5*exp(-x);
outGr = x^4*exp(-x)*(x-5);

end

