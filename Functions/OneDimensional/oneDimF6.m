function [ outVal, outGr ] = oneDimF6( x )
%Function f6 from functions.txt

outVal = sin(x - 1) + x^2*sqrt(abs(x - 1));
outGr = cos(1 - x) + 2*x*sqrt(abs(x - 1)) + (x^2*(x - 1)/2*(abs(x - 1)^(3/2)));

end

