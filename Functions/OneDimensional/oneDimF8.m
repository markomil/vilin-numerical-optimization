function [ outVal, outGr ] = oneDimF8( x )
%Function f8 from functions.txt

outVal = -1 + cos(x)*exp(-sqrt(x));
outGr = -1 * exp(-sqrt(x))*(2*sqrt(x)*sin(x)+cos(x))/2*sqrt(x);

end

