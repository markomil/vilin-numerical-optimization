function [ outVal, outGr ] = oneDimF4( x )
%Function f4 from functions.txt

outVal = -x*exp(-x^2);
outGr = exp(-x^2)*(2*x^2-1);

end

